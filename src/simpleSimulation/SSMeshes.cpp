#include "SSMeshes.h"

#include "simpleSimulation/SSBody.h"

#include <QDebug>
#include <QString>

// ELLIPSOID
casadi::MX SSEllipsoidMesh::constraintJacobian(casadi::MX gamma, casadi::MX q) {
    //qDebug() << "Using manual Jacobi (Ellipsoid)!";
    using namespace casadi;
    // 1. Parameter entpacken (Exakt wie in der Distanzfunktion)
    MX p_center = q(Slice(0, 3));
    MX p_rot = q(Slice(3, 12));
    MX R = MX::reshape(p_rot, 3, 3);
    
    // 2. Relativvektor im lokalen Koordinatensystem
    // x_local(0) entspricht ((gamma - phi) * d1)
    // x_local(1) entspricht ((gamma - phi) * d2)
    // x_local(2) entspricht ((gamma - phi) * d3)
    MX x_local = MX::mtimes(R.T(), (gamma - p_center));

    // 3. Direktoren (Achsen) extrahieren -> Das sind die Spalten von R
    MX d1 = R(Slice(), 0);
    MX d2 = R(Slice(), 1);
    MX d3 = R(Slice(), 2);

    // 4. Ellipsoid-Jacobian nach Gleichung (12):
    MX term1 = (2.0 * x_local(0) / (A * A)) * d1;
    MX term2 = (2.0 * x_local(1) / (B * B)) * d2;
    MX term3 = (2.0 * x_local(2) / (C * C)) * d3;
    
    return term1 + term2 + term3;
    
}
    
casadi::MX SSEllipsoidMesh::constraintDistance(casadi::MX gamma, casadi::MX q) {
    using namespace casadi;
    // 1. Parameter entpacken
    // q: [Pos(3), RotCol1(3), RotCol2(3), RotCol3(3)] -> 12 Elemente
    MX p_center = q(Slice(0, 3));

    /* MX v1 = q(Slice(3, 6));  // Achsenvektor 1 (lokale x-Achse)
    MX v2 = q(Slice(6, 9));  // Achsenvektor 2 (lokale y-Achse)
    MX v3 = q(Slice(9, 12)); // Achsenvektor 3 (lokale z-Achse)  
    MX v1 = MX::vertcat({q(3), q(6), q(9)}); 
    MX v2 = MX::vertcat({q(4), q(7), q(10)}); 
    MX v3 = MX::vertcat({q(5), q(8), q(11)});

    // 2. Relativvektor berechnen (x im alten Code)
    MX x_rel = gamma - p_center;
    // 3. Projektionen auf die Achsen (Dot Product)
    MX proj1 = MX::mtimes(x_rel.T(), v1);
    MX proj2 = MX::mtimes(x_rel.T(), v2);
    MX proj3 = MX::mtimes(x_rel.T(), v3);
    */
    MX p_rot = q(Slice(3, 12)); // Holt die 9 Rotationswerte
    MX R = MX::reshape(p_rot, 3, 3); // R ist jetzt korrekt rekonstruiert
    // Für Distanzberechnung (Global -> Lokal) brauchst du R^T
    MX x_local = MX::mtimes(R.T(), (gamma - p_center));

    // 4. Ellipsoid-Gleichung: (proj/a)^2 + ... - 1
    // Wir nutzen hier die Member-Variablen A, B, C deiner Klasse
    //qDebug() << "A: " << A << ", B: " << B << ", C: " << C;
    /* MX term1 = (proj1 * proj1) / (A * A);
    MX term2 = (proj2 * proj2) / (B * B);
    MX term3 = (proj3 * proj3) / (C * C); */
    MX term1 = (x_local(0) * x_local(0)) / (A * A);
    MX term2 = (x_local(1) * x_local(1)) / (B * B);
    MX term3 = (x_local(2) * x_local(2)) / (C * C);
    MX h = term1 + term2 + term3 - 1.0;
    return h; 
    
}


casadi::MX SSEllipsoidMesh::constraintJacobianLocal(casadi::MX gamma, casadi::MX q)
{
    using namespace casadi;

    // 1. Center Position
    MX p_center = q(Slice(0, 3));

    // 2. Rotations-Matrix Elemente einzeln holen (da p Zeilenweise gespeichert ist)
    // Wir brauchen aber die SPALTEN-Vektoren für die Rücktransformation des Gradienten!
    // Spalte 0 (Lokale X-Achse im Globalen Raum)
    MX col0 = MX::vertcat({q(3), q(6), q(9)}); 
    // Spalte 1 (Lokale Y-Achse im Globalen Raum)
    MX col1 = MX::vertcat({q(4), q(7), q(10)}); 
    // Spalte 2 (Lokale Z-Achse im Globalen Raum)
    MX col2 = MX::vertcat({q(5), q(8), q(11)}); 

    // Für die Projektion (Global -> Lokal) brauchen wir die ZEILEN (oder R transponiert)
    // Zeile 0 (entspricht Spalte 0 von R^T)
    MX row0 = MX::vertcat({q(3), q(4), q(5)});
    MX row1 = MX::vertcat({q(6), q(7), q(8)});
    MX row2 = MX::vertcat({q(9), q(10), q(11)});

    // 3. Relativvektor (Global)
    MX x_rel_global = gamma - p_center;

    // 4. Projektion in den lokalen Raum (x_local = R^T * x_global)
    // Dazu projizieren wir auf die ZEILEN von R
    MX x_loc_x = MX::mtimes(x_rel_global.T(), row0);
    MX x_loc_y = MX::mtimes(x_rel_global.T(), row1);
    MX x_loc_z = MX::mtimes(x_rel_global.T(), row2);

    // 5. Gradient im lokalen Raum berechnen ( dF/dx_loc = 2*x_loc / axis^2 )
    MX grad_loc_x = (2.0 * x_loc_x / (A * A));
    MX grad_loc_y = (2.0 * x_loc_y / (B * B));
    MX grad_loc_z = (2.0 * x_loc_z / (C * C));

    // 6. Gradient in den globalen Raum zurücktransformieren ( dF/dx_glob = R * dF/dx_loc )
    // Dazu multiplizieren wir die lokalen Anteile mit den SPALTEN von R (den Achsen)
    MX grad_global = grad_loc_x * col0 + grad_loc_y * col1 + grad_loc_z * col2;

    return grad_global;
}

casadi::MX SSEllipsoidMesh::constraintDistanceLocal(casadi::MX gamma, casadi::MX q)
{
    using namespace casadi;

    // 1. Parameter entpacken
    // q: [Pos(3), RotCol1(3), RotCol2(3), RotCol3(3)] -> 12 Elemente
    MX p_center = q(Slice(0, 3));
    MX v1 = q(Slice(3, 6));  // Achsenvektor 1 (lokale x-Achse)
    MX v2 = q(Slice(6, 9));  // Achsenvektor 2 (lokale y-Achse)
    MX v3 = q(Slice(9, 12)); // Achsenvektor 3 (lokale z-Achse)
    //MX v1 = MX::vertcat({q(3), q(6), q(9)}); // Spalte 0 (Lokale X-Achse): Indizes 3, 6, 9 -> (R00, R10, R20)
    //MX v2 = MX::vertcat({q(4), q(7), q(10)}); // Spalte 1 (Lokale Y-Achse): Indizes 4, 7, 10 -> (R01, R11, R21)
    //MX v3 = MX::vertcat({q(5), q(8), q(11)}); // Spalte 2 (Lokale Z-Achse): Indizes 5, 8, 11 -> (R02, R12, R22)
    // 2. Relativvektor berechnen (x im alten Code)
    MX x_rel = gamma - p_center;
    // 3. Projektionen auf die Achsen (Dot Product)
    MX proj1 = MX::mtimes(x_rel.T(), v1);
    MX proj2 = MX::mtimes(x_rel.T(), v2);
    MX proj3 = MX::mtimes(x_rel.T(), v3);
    // 4. Ellipsoid-Gleichung: (proj/a)^2 + ... - 1
    // Wir nutzen hier die Member-Variablen A, B, C deiner Klasse
    MX term1 = (proj1 * proj1) / (A * A);
    MX term2 = (proj2 * proj2) / (B * B);
    MX term3 = (proj3 * proj3) / (C * C);
    MX h = term1 + term2 + term3 - 1.0;
    return h; 
}

double SSEllipsoidMesh::getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance)
{
    if (signedDistance == false) { 
        // Einfache euklidische Distanz zum Zentrum
        if (GlobalDiscreteMeshPoints.empty()) {
            // 1. In lokales Koordinatensystem des Ellipsoids transformieren
            MWMath::Point3D diff = pGlobal - PositionGlobal;
            MWMath::Point3D pLoc = OrientationGlobal.transposed().transform(diff);
            // Fallback: Wenn keine diskreten Punkte vorhanden sind, verwenden wir den Abstand zum Zentrum
            return std::sqrt(pLoc.x * pLoc.x + pLoc.y * pLoc.y + pLoc.z * pLoc.z);
        }
        else{
            //qDebug() << "    -> Using discrete mesh points for distance calculation.";
            // 2. Nächsten Punkt im diskreten Mesh finden
            double minDistSq = std::numeric_limits<double>::max();
            for (const auto& discMeshP : GlobalDiscreteMeshPoints) {
                double distSq = MWMath::distance(pGlobal, discMeshP);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                }
            }
            return minDistSq;
        }
    } 
    else {
        // 1. In lokales Koordinatensystem des Ellipsoids transformieren
        MWMath::Point3D diff = pGlobal - PositionGlobal;
        MWMath::Point3D pLoc = OrientationGlobal.transposed().transform(diff);
        // Algebraische SDF für Ellipsoid: (x/a)² + (y/b)² + (z/c)² - 1
        // Wert < 0: Punkt ist innen (Kollision)
        // Wert > 0: Punkt ist außen
        double val = (pLoc.x * pLoc.x) / (A * A) + 
                     (pLoc.y * pLoc.y) / (B * B) + 
                     (pLoc.z * pLoc.z) / (C * C);
        
        return val - 1.0; 
    }
}

void SSEllipsoidMesh::discretizeMesh(int discrCount)
{
    // 1. Alten Inhalt löschen, um sauberen Start zu haben
    int expectedSize = (discrCount+1) * discrCount;
    if(GlobalDiscreteMeshPoints.empty() || GlobalDiscreteMeshPoints.size() != expectedSize){
        GlobalDiscreteMeshPoints.resize(expectedSize);
    }

    // Sicherheitscheck: Zu geringe Auflösung vermeiden
    if (discrCount < 3) discrCount = 3;

    for (int i = 0; i <= discrCount; ++i) {
        
        // normalizedPos geht von -1.0 bis +1.0
        double normalizedPos = -1.0 + 2.0 * (double)i / (double)discrCount;
        
        // Lokale Z-Koordinate (Höhe)
        double zLoc = normalizedPos * C;

        // Bei z=0 ist der Faktor 1 (voller Äquator), bei z=C ist er 0 (Pol).
        double term = 1.0 - (normalizedPos * normalizedPos);
        
        // Wurzel aus negativer Zahl verhindern
        if (term < 0.0) term = 0.0;
        
        double radiusFactor = std::sqrt(term);

        // --- SCHRITT C: Um die Z-Achse rotieren (Ring erzeugen) ---
        // 360° in discrCount Schritte aufteilen
        for (int j = 0; j < discrCount; ++j) {
            double angle = (2.0 * M_PI * (double)j) / (double)discrCount;

            double xLoc = A * radiusFactor * std::cos(angle);
            double yLoc = B * radiusFactor * std::sin(angle);

            // lokalen Koordinatensystem 
            MWMath::Point3D pLoc(xLoc, yLoc, zLoc);

            MWMath::Point3D pGlob = PositionGlobal + OrientationGlobal.transform(pLoc);
            GlobalDiscreteMeshPoints[i * discrCount + j] = pGlob;
        }
    }
}

// CYLINDER
casadi::MX SSCylinderMesh::constraintDistance(casadi::MX gamma, casadi::MX q) {
    casadi::MX phi = q(casadi::Slice(0, 3));
    casadi::MX d1  = q(casadi::Slice(3, 6));
    casadi::MX d2  = q(casadi::Slice(6, 9));
    casadi::MX diff = gamma - phi;

    // dist = ((diff' * d1)^2 / R^2) + ((diff' * d2)^2 / R^2) - 1
    return (pow(casadi::MX::dot(diff, d1), 2) / (Radius * Radius)) + 
           (pow(casadi::MX::dot(diff, d2), 2) / (Radius * Radius)) - 1.0;
}

double SSCylinderMesh::getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance) {
    // 1. In lokales System transformieren
    MWMath::Point3D diff = pGlobal - PositionGlobal;
    MWMath::Point3D pLoc = OrientationGlobal.transposed().transform(diff); 

    // Wir nehmen an: d1, d2 spannen Radius auf, d3 ist Längsachse (Z)
    double r_xy = std::sqrt(pLoc.x * pLoc.x + pLoc.y * pLoc.y);
    double halfHeight = Height / 2.0;

    // Distanzberechnung (SDF - Signed Distance Function Logik)
    double dist_r = r_xy - Radius;
    double dist_z = std::abs(pLoc.z) - halfHeight;

    if (dist_z > 0 && dist_r > 0) {
        // Punkt ist diagonal außerhalb der Kante
        return std::sqrt(dist_r * dist_r + dist_z * dist_z);
    }
    // Sonst das Maximum der Abweichungen (Standard Zylinder-Mantel oder flache Kappe)
    return std::max(dist_r, dist_z);
}

casadi::MX SSCylinderMesh::constraintJacobian(casadi::MX gamma, casadi::MX q) {
    //qDebug() << "Using manual Jacobi (Cylinder)!";
    using namespace casadi;
    
    // 1. Parameter entpacken (Exakt wie in der Distanzfunktion)
    MX p_center = q(Slice(0, 3));
    MX p_rot = q(Slice(3, 12));
    MX R = MX::reshape(p_rot, 3, 3);
    
    // 2. Relativvektor im lokalen Koordinatensystem
    MX x_local = MX::mtimes(R.T(), (gamma - p_center));

    // 3. Direktoren (Achsen) extrahieren -> Wir brauchen für den Zylinder nur d1 und d2
    MX d1 = R(Slice(), 0);
    MX d2 = R(Slice(), 1);

    // 4. Cylinder-Jacobian nach Gleichung (14):
    // ACHTUNG: Keine Teilung durch den Radius hier!
    MX term1 = 2.0 * x_local(0) * d1;
    MX term2 = 2.0 * x_local(1) * d2;
    
    return term1 + term2;
}


// TORUS
casadi::MX SSTorusMesh::constraintDistance(casadi::MX gamma, casadi::MX q) {
    casadi::MX phi = q(casadi::Slice(0, 3));
    casadi::MX d1  = q(casadi::Slice(3, 6));
    casadi::MX d2  = q(casadi::Slice(6, 9));
    casadi::MX d3  = q(casadi::Slice(9, 12));
    casadi::MX diff = gamma - phi;

    // Skalierte Projektionen (Inverse Skalierung)
    casadi::MX dot1 = casadi::MX::dot(diff, d1) / A;
    casadi::MX dot2 = casadi::MX::dot(diff, d2) / B;
    casadi::MX dot3 = casadi::MX::dot(diff, d3) / C;
    /* casadi::MX dot1 = casadi::MX::dot(diff, d1);
    casadi::MX dot2 = casadi::MX::dot(diff, d2);
    casadi::MX dot3 = casadi::MX::dot(diff, d3); */

    // Wir nutzen die Einheits-Torus-Gleichung (R und r beziehen sich auf den skalierten Raum)
    // dist = (dot1^2 + dot2^2 + dot3^2 + R^2 - r^2)^2 - 4*R^2*(dot1^2 + dot2^2)
    casadi::MX term1 = pow(dot1, 2) + pow(dot2, 2) + pow(dot3, 2) + pow(R, 2) - pow(r, 2);
    return pow(term1, 2) - 4.0 * pow(R, 2) * (pow(dot1, 2) + pow(dot2, 2));
}

casadi::MX SSTorusMesh::constraintJacobian(casadi::MX gamma, casadi::MX q) {
    //qDebug() << "Using manual Jacobi (Torus)";
    using namespace casadi;
    
    // 1. Parameter entpacken
    MX p_center = q(Slice(0, 3));
    MX p_rot = q(Slice(3, 12));
    MX R_mat = MX::reshape(p_rot, 3, 3);
    
    // 2. Relativvektor im lokalen Koordinatensystem
    // x_local(0) = lokales X
    // x_local(1) = lokales Y
    // x_local(2) = lokales Z
    MX x_local = MX::mtimes(R_mat.T(), (gamma - p_center));

    // 3. Direktoren (Achsen) extrahieren -> Spalten von R_mat
    MX d1 = R_mat(Slice(), 0);
    MX d2 = R_mat(Slice(), 1);
    MX d3 = R_mat(Slice(), 2);

    // 4. Torus-Jacobian berechnen
    // Abstand des Punktes zur Z-Achse in der X-Y-Ebene: S = sqrt(x^2 + y^2)
    // (+ 1e-12 verhindert Division durch 0, falls der Punkt exakt auf der Z-Achse liegt)
    MX S = MX::sqrt(x_local(0) * x_local(0) + x_local(1) * x_local(1) + 1e-12);
    
    // Vorfaktor für die X- und Y-Ableitung: 2 * (1 - R_major / S)
    // Wir nehmen an, dass R in deiner Klasse der Major-Radius ist (Ring-Radius)
    MX factor_xy = 2.0 * (1.0 - R / S);

    // Ableitungen in die jeweiligen Richtungen
    MX term1 = (factor_xy * x_local(0)) * d1;
    MX term2 = (factor_xy * x_local(1)) * d2;
    MX term3 = (2.0 * x_local(2)) * d3;
    
    return term1 + term2 + term3;
}

double SSTorusMesh::getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance) {
    /* // 1. In lokales System transformieren
    MWMath::Point3D diff = pGlobal - PositionGlobal;
    MWMath::Point3D pLoc = OrientationGlobal.transposed().transform(diff);

    // 2. Punkt in den Einheitsraum skalieren
    double px = pLoc.x / A;
    double py = pLoc.y / B;
    double pz = pLoc.z / C;

    // 3. Distanz im skalierten Raum berechnen
    double dist_xy = std::sqrt(px * px + py * py);
    double delta_R = dist_xy - R;
    double dist_to_center_line = std::sqrt(delta_R * delta_R + pz * pz);
    double d_model = dist_to_center_line - r;

    // 4. Zurückskalieren der Distanz (Näherung durch minimalen Skalierungsfaktor)
    // Das verhindert, dass der Solver "denkt", er sei schon vorbei, obwohl er noch im Mesh ist.
    double min_scale = std::min({A, B, C});
    return d_model * min_scale; */
    if (GlobalDiscreteMeshPoints.empty()) {
        MWMath::Point3D diff = pGlobal - PositionGlobal;
        MWMath::Point3D pLoc = OrientationGlobal.transposed().transform(diff);
        // Fallback: Wenn keine diskreten Punkte vorhanden sind, verwenden wir den Abstand zum Zentrum
        return std::sqrt(pLoc.x * pLoc.x + pLoc.y * pLoc.y + pLoc.z * pLoc.z);
    }
    else{
        double minDistSq = std::numeric_limits<double>::max();
            for (const auto& discMeshP : GlobalDiscreteMeshPoints) {
                double distSq = MWMath::distance(pGlobal, discMeshP);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                }
            }
            return minDistSq;
    }
}

void SSTorusMesh::discretizeMesh(int discrCount)
{
    // Sicherheitscheck
    if (discrCount < 4) discrCount = 4;

    // 1. Größe berechnen
    // Wir laufen 2x im Kreis (Ring und Rohr), daher discrCount * discrCount Punkte
    int expectedSize = discrCount * discrCount;
    
    if(GlobalDiscreteMeshPoints.empty() || GlobalDiscreteMeshPoints.size() != expectedSize){
        GlobalDiscreteMeshPoints.resize(expectedSize);
    }

    // 2. Schleifen über die Winkel
    for (int i = 0; i < discrCount; ++i) {
        
        // Winkel u: "Longitude" (Läuft um den großen Ring / Z-Achse)
        double u = (2.0 * M_PI * (double)i) / (double)discrCount;
        double cos_u = std::cos(u);
        double sin_u = std::sin(u);

        for (int j = 0; j < discrCount; ++j) {
            
            // Winkel v: "Latitude" (Läuft um den Querschnitt des Rohrs)
            double v = (2.0 * M_PI * (double)j) / (double)discrCount;
            
            // --- Parametrisierung des Torus (liegend in XY-Ebene) ---
            // x = (R + r*cos(v)) * cos(u)
            // y = (R + r*cos(v)) * sin(u)
            // z = r * sin(v)
            
            // Abstand vom Torus-Zentrum zum Punkt im Rohr-Querschnitt (in der XY-Ebene)
            double crossSectionRadius = R + r * std::cos(v);

            double xLoc = crossSectionRadius * cos_u;
            double yLoc = crossSectionRadius * sin_u;
            double zLoc = r * std::sin(v);

            // Lokaler Punkt
            MWMath::Point3D pLoc(xLoc, yLoc, zLoc);

            // Transformation in Weltkoordinaten
            // WICHTIG: OrientationGlobal sorgt dafür, dass der Torus im Raum richtig gedreht ist
            MWMath::Point3D pGlob = PositionGlobal + OrientationGlobal.transform(pLoc);
            
            // Speichern (linearisiertes 2D Array)
            GlobalDiscreteMeshPoints[i * discrCount + j] = pGlob;
        }
    }
}

// SSMESH
void SSMesh::InitializeMesh()
{
    if (Parent){
        Parent->Meshes.push_back(shared_from_this());
        if (Parent->bDebug) qDebug() << "Mesh " << QString::fromStdString(Name) << " added to Parent " << QString::fromStdString(Parent->Name);
        
        if (dynamic_cast<SSJoint*>(Parent.get()))// != nullptr && dynamic_cast<SSTorusMesh*>(shared_from_this().get()) == nullptr)
        {
            //qDebug() << "  -> Parent is a Joint.";
            bIsJointMesh = true;
        }
    } // Füge dieses Mesh der Liste des Parents hinzu

}

void SSMesh::getCasadiParentGPOs(casadi::MX &pos, casadi::MX &ori)
{
        MWMath::Point3D positionGlobal;
        MWMath::RotMatrix3x3 orientationGlobal;
        if (Parent) {
            positionGlobal = Parent->PositionGlobal;
            orientationGlobal = Parent->OrientationGlobal;
        }
        else {
            positionGlobal = Position2ParentRelInParentFrame;
            orientationGlobal = Orientation2ParentRel;
        }
        pos = casadi::MX::vertcat({positionGlobal.x, positionGlobal.y, positionGlobal.z});
        ori = casadi::MX::zeros(3,3);
        for (int i=0;i<3;i++){for (int j=0;j<3;j++){ori(i,j) = orientationGlobal.m[i][j];}}
}

void SSMesh::updateMeshPosAndRot() {
    // Aktualisiere globale Position und Rotation basierend auf Parent
    // qDebug() << "Updating Mesh Position: " <<  QString::fromStdString(Name);
    if (Parent) {
        // qDebug() << "  old: Pos: " << QString::fromStdString(PositionGlobal.print()) << ", Ori: " << QString::fromStdString(OrientationGlobal.print());
        PositionGlobal = Parent->PositionGlobal + Parent->OrientationGlobal.transform(Position2ParentRelInParentFrame);
        OrientationGlobal = Parent->OrientationGlobal * Orientation2ParentRel;
        //PositionGlobal = Parent->PositionGlobal + OrientationGlobal.transform(Position2ParentRelInParentFrame);
        if (Name == "_"){
            qDebug() << "    " << this->Name.c_str() << ": Pos: " << QString::fromStdString(PositionGlobal.print()) << ", Ori: " << QString::fromStdString(OrientationGlobal.print());
        }
    } else {
        PositionGlobal = Position2ParentRelInParentFrame; // Falls kein Parent, dann ist die lokale Position die globale
        OrientationGlobal = Orientation2ParentRel; // Falls kein Parent, dann ist die lokale Rotation die globale
    }
    qDebug() << "Global Position: " << Name.c_str() << "=" <<  PositionGlobal.print().c_str();

    
};


// ELLIPTICAL TORUS
casadi::MX SSEllipticalTorusMesh::constraintJacobian(casadi::MX gamma, casadi::MX q)
{
    casadi::MX phi = q(casadi::Slice(0, 3));
    casadi::MX d1  = q(casadi::Slice(3, 6));
    casadi::MX d2  = q(casadi::Slice(6, 9));
    casadi::MX d3  = q(casadi::Slice(9, 12));

    casadi::MX diff = gamma - phi;

    // Lokale Koordinaten
    casadi::MX u = casadi::MX::dot(diff, d1);
    casadi::MX v = casadi::MX::dot(diff, d2);
    casadi::MX w = casadi::MX::dot(diff, d3);

    // Hilfsterme (identisch zur Distance)
    casadi::MX term_path = (u * u) / (RX * RX) + (v * v) / (RY * RY);
    casadi::MX sqrt_path = sqrt(term_path + 1e-16);

    // --- Partielle Ableitungen (Kettenregel) ---
    
    // Wir leiten den ersten Term T1 nach u und v ab.
    // T1 = ( (sqrt(S) - 1) * RX / R1 )^2
    // dT1/dS = 2 * ( ... ) * (RX / R1) * (1 / (2 * sqrt(S)))
    //        = ( (sqrt(S) - 1) * RX / R1 ) * (RX / R1) * (1 / sqrt(S))
    
    // Klammer-Ausdruck (Radialer Fehler normiert auf r1)
    casadi::MX pre_factor = ((sqrt_path - 1.0) * RX / R1);
    
    // Zusammenfassung der Faktoren für die Ableitung nach S
    casadi::MX d_term1_dS = 2.0 * pre_factor * (RX / R1) * (0.5 / sqrt_path);
    // Vereinfacht:
    // casadi::MX d_term1_dS = pre_factor * (RX / R1) / sqrt_path;

    // Ableitung von S nach u und v
    casadi::MX dS_du = 2.0 * u / (RX * RX);
    casadi::MX dS_dv = 2.0 * v / (RY * RY);

    // Totale Ableitung nach u und v
    casadi::MX df_du = d_term1_dS * dS_du;
    casadi::MX df_dv = d_term1_dS * dS_dv;

    // Ableitung des zweiten Terms T2 = (w / R2)^2 nach w
    casadi::MX df_dw = 2.0 * w / (R2 * R2);

    // --- Rücktransformation ---
    // Gradient = df_du * d1 + df_dv * d2 + df_dw * d3
    casadi::MX grad = df_du * d1 + df_dv * d2 + df_dw * d3;

    return grad;
}

casadi::MX SSEllipticalTorusMesh::constraintDistance(casadi::MX gamma, casadi::MX q)
{
    // Parameter extrahieren
    casadi::MX phi = q(casadi::Slice(0, 3));
    casadi::MX d1  = q(casadi::Slice(3, 6));
    casadi::MX d2  = q(casadi::Slice(6, 9));
    casadi::MX d3  = q(casadi::Slice(9, 12));

    casadi::MX diff = gamma - phi;

    // --- 1. Projektion ins lokale System ---
    casadi::MX u = casadi::MX::dot(diff, d1);
    casadi::MX v = casadi::MX::dot(diff, d2);
    casadi::MX w = casadi::MX::dot(diff, d3);

    // --- 2. Berechnung des Pfad-Terms (Große Ellipse) ---
    // Wir berechnen, wo wir relativ zur Pfad-Ellipse stehen.
    // term_path = 1.0 bedeutet: Wir sind exakt auf der Pfad-Linie.
    casadi::MX term_path = (u * u) / (RX * RX) + (v * v) / (RY * RY);
    casadi::MX sqrt_path = sqrt(term_path + 1e-16); // epsilon für Stabilität

    // Der Term (sqrt_path - 1.0) ist der "Abstand" zur Pfad-Ellipse.
    // Aber Achtung: Er ist dimensionslos!
    // Um ihn mit r1 (Meter) zu vergleichen, skalieren wir ihn mit RX.
    // (Das ist eine Approximation; exakt wäre es nur bei einem Kreis).
    casadi::MX radial_dist_approx = (sqrt_path - 1.0) * RX;

    // --- 3. Zusammensetzen mit Querschnitts-Ellipse (Kleine Ellipse) ---
    // Formel: (dist_radial / r1)^2 + (dist_vertical / r2)^2 - 1
    
    casadi::MX term_r1 = (radial_dist_approx * radial_dist_approx) / (R1 * R1);
    casadi::MX term_r2 = (w * w) / (R2 * R2);

    return term_r1 + term_r2 - 1.0;
}

double SSEllipticalTorusMesh::getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance)
{
    if (GlobalDiscreteMeshPoints.empty()) {
        MWMath::Point3D diff = pGlobal - PositionGlobal;
        MWMath::Point3D pLoc = OrientationGlobal.transposed().transform(diff);
        // Fallback: Wenn keine diskreten Punkte vorhanden sind, verwenden wir den Abstand zum Zentrum
        return std::sqrt(pLoc.x * pLoc.x + pLoc.y * pLoc.y + pLoc.z * pLoc.z);
    }
    else{
        double minDistSq = std::numeric_limits<double>::max();
            for (const auto& discMeshP : GlobalDiscreteMeshPoints) {
                double distSq = MWMath::distance(pGlobal, discMeshP);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                }
            }
            return minDistSq;
    }
}

void SSEllipticalTorusMesh::discretizeMesh(int discrCount)
{
    discrCount = 31; // TEST
    // 1. Sicherheitscheck
    if (discrCount < 4) discrCount = 4;

    // Größe berechnen (discCount * discCount Punkte)
    int expectedSize = discrCount * discrCount;
    
    if(GlobalDiscreteMeshPoints.size() != expectedSize){
        GlobalDiscreteMeshPoints.resize(expectedSize);
    }

    // 2. Schleifen über die Winkel
    for (int i = 0; i < discrCount; ++i) {
        
        // Winkel u: Parameter der Pfad-Ellipse (um die Z-Achse)
        double u = (2.0 * M_PI * (double)i) / (double)discrCount;
        double cos_u = std::cos(u);
        double sin_u = std::sin(u);

        // Punkt auf der Pfad-Ellipse (Mittellinie des Schlauchs)
        // x_path = RX * cos(u)
        // y_path = RY * sin(u)

        for (int j = 0; j < discrCount; ++j) {
            
            // Winkel v: Parameter der Rohr-Ellipse (Querschnitt)
            double v = (2.0 * M_PI * (double)j) / (double)discrCount;
            double cos_v = std::cos(v);
            double sin_v = std::sin(v);

            // --- Parametrisierung des elliptischen Torus ---
            // Wir berechnen zuerst den Punkt im lokalen Koordinatensystem:
            
            // Der Radius der Rohr-Ellipse wirkt sich auf die Verschiebung aus.
            // Der Querschnitt "schwillt" um den Pfad herum an.
            // xLoc und yLoc nutzen die Pfad-Radien RX/RY + Rohr-Anteil
            // zLoc nutzt den vertikalen Rohr-Radius R2
            
            double xLoc = (RX + R1 * cos_v) * cos_u;
            double yLoc = (RY + R1 * cos_v) * sin_u;
            double zLoc = R2 * sin_v;

            // Lokaler Punkt
            MWMath::Point3D pLoc(xLoc, yLoc, zLoc);

            // 3. Transformation in Weltkoordinaten
            // Nutzt die von SSBody/SSMesh berechnete globale Pose
            MWMath::Point3D pGlob = PositionGlobal + OrientationGlobal.transform(pLoc);
            
            // Speichern
            GlobalDiscreteMeshPoints[i * discrCount + j] = pGlob;
        }
    }
}


// SMILE MESH
casadi::MX SSSmileMesh::constraintDistance(casadi::MX gamma, casadi::MX q) {
    casadi::MX phi  = q(casadi::Slice(0, 3));
    casadi::MX d1   = q(casadi::Slice(3, 6));
    casadi::MX d2   = q(casadi::Slice(6, 9));
    casadi::MX d3   = q(casadi::Slice(9, 12));
    casadi::MX diff = gamma - phi;

    // Skalierte Projektionen in den lokalen Raum
    casadi::MX x = casadi::MX::dot(diff, d1) / A;
    casadi::MX y = casadi::MX::dot(diff, d2) / B;
    casadi::MX z = casadi::MX::dot(diff, d3) / C;

    // Quadrate vorberechnen, da sie mehrfach gebraucht werden
    casadi::MX x2 = pow(x, 2);
    casadi::MX y2 = pow(y, 2);
    casadi::MX z2 = pow(z, 2);

    // Basis-Terme der Smile-Surface
    casadi::MX U = y - x2 - y2 + 1.0;
    casadi::MX V = x2 + y2 + z2;

    // Distanzfunktion (impliziter Wert)
    return pow(U, 4) + pow(V, 4) - 1.0;
}

void SSSmileMesh::discretizeMesh(int discrCount)
{
    // Sicherheitscheck
    if (discrCount < 4) discrCount = 4;

    // Da wir vorab nicht wissen, wie viele Punkte exakt gültig sind (Wurzel-Bedingung) 
    // und wir zudem eine obere und untere Hälfte haben, nutzen wir einen dynamischen Aufbau.
    GlobalDiscreteMeshPoints.clear();
    
    // Wir reservieren grob Speicher, um Reallokationen zu vermeiden (ca. 2x für +Z und -Z)
    GlobalDiscreteMeshPoints.reserve(discrCount * discrCount * 2);

    // 1. Schleife: Winkel theta (in der XY-Ebene, von 0 bis 2*PI)
    for (int i = 0; i < discrCount; ++i) {
        double theta = (2.0 * M_PI * (double)i) / (double)discrCount;
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);

        // 2. Schleife: Parameter w (von -1 bis 1)
        for (int j = 0; j < discrCount; ++j) {
            
            // Wir mappen j von [0, discrCount-1] auf w in [-1.0, 1.0]
            double w = -1.0 + (2.0 * (double)j) / (double)(discrCount - 1);
            
            // --- Schritt A: Radius rho berechnen ---
            // Formel: rho = (sin(theta) + sqrt(sin^2(theta) + 4(1 - w))) / 2
            double term_under_rho_root = (sin_theta * sin_theta) + 4.0 * (1.0 - w);
            
            // Sicherheitshalber abfangen (sollte wegen w <= 1 immer positiv sein)
            if (term_under_rho_root < 0.0) continue; 
            
            double rho = (sin_theta + std::sqrt(term_under_rho_root)) / 2.0;

            // --- Schritt B: Z-Koordinate berechnen ---
            // Bedingung: z^2 = sqrt[4](1 - w^4) - rho^2
            double w4 = w * w * w * w;
            double B_term = std::pow(1.0 - w4, 0.25); // 4. Wurzel aus (1 - w^4)
            double rho2 = rho * rho;

            double term_under_z_root = B_term - rho2;

            // Nur wenn der Term positiv ist, gehört der Punkt (w, theta) zur echten Oberfläche
            if (term_under_z_root >= 0.0) {
                
                // Basis-Koordinaten unskaliert
                double z_base = std::sqrt(term_under_z_root);
                double x_base = rho * cos_theta;
                double y_base = rho * sin_theta;

                // Mit A, B, C skalieren
                double xLoc = A * x_base;
                double yLoc = B * y_base;
                double zLoc_pos = C * z_base;
                double zLoc_neg = C * -z_base;

                // --- Obere Schale (+Z) ---
                MWMath::Point3D pLoc_pos(xLoc, yLoc, zLoc_pos);
                MWMath::Point3D pGlob_pos = PositionGlobal + OrientationGlobal.transform(pLoc_pos);
                GlobalDiscreteMeshPoints.push_back(pGlob_pos);

                // --- Untere Schale (-Z) ---
                // Wir fügen die negative Seite nur hinzu, wenn z nicht exakt 0 ist, 
                // um doppelte Punkte an den Nahtstellen/Kanten zu vermeiden.
                if (z_base > 1e-6) {
                    MWMath::Point3D pLoc_neg(xLoc, yLoc, zLoc_neg);
                    MWMath::Point3D pGlob_neg = PositionGlobal + OrientationGlobal.transform(pLoc_neg);
                    GlobalDiscreteMeshPoints.push_back(pGlob_neg);
                }
            }
        }
    }
    
    // Optional: Überflüssigen Speicher freigeben
    GlobalDiscreteMeshPoints.shrink_to_fit();
}

double SSSmileMesh::getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance) {
    if (GlobalDiscreteMeshPoints.empty()) {
        MWMath::Point3D diff = pGlobal - PositionGlobal;
        MWMath::Point3D pLoc = OrientationGlobal.transposed().transform(diff);
        // Fallback: Wenn keine diskreten Punkte vorhanden sind, verwenden wir den Abstand zum Zentrum.
        // Bei skalierten Flächen ist das sehr grob, aber als Fallback absolut ausreichend.
        return std::sqrt(pLoc.x * pLoc.x + pLoc.y * pLoc.y + pLoc.z * pLoc.z);
    }
    else {
        double minDistSq = std::numeric_limits<double>::max();
        for (const auto& discMeshP : GlobalDiscreteMeshPoints) {
            // Berechnet die Distanz (oder quadrierte Distanz, je nach Implementierung deiner Library)
            double dist = MWMath::distance(pGlobal, discMeshP);
            // Wir gehen hier davon aus, dass wir Vergleiche optimieren (falls distSq genutzt wird)
            if (dist < minDistSq) {
                minDistSq = dist;
            }
        }
        
        // Falls deine Funktion MWMath::distance tatsächlich die quadrierte Distanz zurückgibt,
        // solltest du hier return std::sqrt(minDistSq); verwenden. 
        // Wenn sie die echte Distanz liefert, reicht return minDistSq;
        return minDistSq; 
    }
}

casadi::MX SSSmileMesh::constraintJacobian(casadi::MX gamma, casadi::MX q) {
    //qDebug() << "Using manual Jacobi (Smile Surface)";
    using namespace casadi;
    
    // 1. Parameter entpacken
    MX p_center = q(Slice(0, 3));
    MX p_rot    = q(Slice(3, 12));
    MX R_mat    = MX::reshape(p_rot, 3, 3);
    
    // 2. Relativvektor zum Zentrum
    MX diff = gamma - p_center;

    // 3. Direktoren (Achsen) extrahieren -> Spalten von R_mat
    MX d1 = R_mat(Slice(), 0);
    MX d2 = R_mat(Slice(), 1);
    MX d3 = R_mat(Slice(), 2);

    // Skalierte lokale Koordinaten ermitteln
    MX x = MX::dot(diff, d1) / A;
    MX y = MX::dot(diff, d2) / B;
    MX z = MX::dot(diff, d3) / C;

    MX x2 = pow(x, 2);
    MX y2 = pow(y, 2);
    MX z2 = pow(z, 2);

    // 4. Hilfsterme berechnen (U und V wie in der Distanzfunktion)
    MX U = y - x2 - y2 + 1.0;
    MX V = x2 + y2 + z2;

    MX U3 = pow(U, 3);
    MX V3 = pow(V, 3);

    // 5. Partielle Ableitungen nach den skalierten lokalen Achsen
    // dF/dx = 8x * (V^3 - U^3)
    MX dF_dx = 8.0 * x * (V3 - U3);
    
    // dF/dy = 4U^3 * (1 - 2y) + 8y * V^3
    MX dF_dy = 4.0 * U3 * (1.0 - 2.0 * y) + 8.0 * y * V3;
    
    // dF/dz = 8z * V^3
    MX dF_dz = 8.0 * z * V3;

    // 6. Kettenregel für den Gradienten im globalen Raum
    // Multiplikation mit den Basisvektoren und der inneren Ableitung der Skalierung (1/A, 1/B, 1/C)
    MX term1 = (dF_dx / A) * d1;
    MX term2 = (dF_dy / B) * d2;
    MX term3 = (dF_dz / C) * d3;
    
    return term1 + term2 + term3;
}

