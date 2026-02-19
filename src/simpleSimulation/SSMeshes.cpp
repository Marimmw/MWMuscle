#include "SSMeshes.h"

#include "simpleSimulation/SSBody.h"

#include <QDebug>
#include <QString>

// ELLIPSOID
casadi::MX SSEllipsoidMesh::constraintJacobian(casadi::MX gamma, casadi::MX q) {
    /* using namespace casadi;
    MX p_center = q(Slice(0, 3));
    //MX v1 = q(Slice(3, 6));
    //MX v2 = q(Slice(6, 9));
    //MX v3 = q(Slice(9, 12)); 
    MX v1 = MX::vertcat({q(3), q(6), q(9)}); 
    MX v2 = MX::vertcat({q(4), q(7), q(10)}); 
    MX v3 = MX::vertcat({q(5), q(8), q(11)});
    MX x_rel = gamma - p_center;
    MX proj1 = MX::mtimes(x_rel.T(), v1);
    MX proj2 = MX::mtimes(x_rel.T(), v2);
    MX proj3 = MX::mtimes(x_rel.T(), v3);
    MX g1 = (2.0 * proj1 / (A * A)) * v1;
    MX g2 = (2.0 * proj2 / (B * B)) * v2;
    MX g3 = (2.0 * proj3 / (C * C)) * v3;
    MX grad = g1 + g2 + g3;
    return grad; */
    
    using namespace casadi;
    MX p_center = q(Slice(0, 3));
    MX p_rot    = q(Slice(3, 12));
    // 2. Matrix wiederherstellen
    // R_transposed enthält die lokalen Achsen als ZEILEN (wenn in C++ Zeilenweise gefüllt wurde).
    // Das ist mathematisch R^T.
    MX R_transposed = MX::reshape(p_rot, 3, 3);
    MX x_rel = gamma - p_center;
    MX x_loc = MX::mtimes(R_transposed, x_rel);
    MX g_loc_x = (2.0 * x_loc(0)) / (A * A);
    MX g_loc_y = (2.0 * x_loc(1)) / (B * B);
    MX g_loc_z = (2.0 * x_loc(2)) / (C * C); 
    MX grad_loc = MX::vertcat({g_loc_x, g_loc_y, g_loc_z});
    MX grad_global = MX::mtimes(R_transposed.T(), grad_loc);
    return grad_global;
    
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
    casadi::MX phi = q(casadi::Slice(0, 3));
    casadi::MX d1  = q(casadi::Slice(3, 6));
    casadi::MX d2  = q(casadi::Slice(6, 9));
    casadi::MX diff = gamma - phi;

    // grad = (2/R^2) * ((diff' * d1)*d1 + (diff' * d2)*d2)
    return (2.0 / (Radius * Radius)) * (
        casadi::MX::dot(diff, d1) * d1 + 
        casadi::MX::dot(diff, d2) * d2
    );
}

/*
casadi::MX SSTorusMesh::constraintDistance(casadi::MX gamma, casadi::MX q) {
    casadi::MX phi = q(casadi::Slice(0, 3));
    casadi::MX d1  = q(casadi::Slice(3, 6));
    casadi::MX d2  = q(casadi::Slice(6, 9));
    casadi::MX d3  = q(casadi::Slice(9, 12));
    casadi::MX diff = gamma - phi;

    casadi::MX dot1 = casadi::MX::dot(diff, d1);
    casadi::MX dot2 = casadi::MX::dot(diff, d2);
    casadi::MX dot3 = casadi::MX::dot(diff, d3);

    // dist = (dot1^2 + dot2^2 + dot3^2 + a^2 - b^2)^2 - 4*a^2*(dot1^2 + dot2^2)
    casadi::MX term1 = pow(dot1, 2) + pow(dot2, 2) + pow(dot3, 2) + pow(R, 2) - pow(r, 2);
    return pow(term1, 2) - 4.0 * pow(R, 2) * (pow(dot1, 2) + pow(dot2, 2));
}

double SSTorusMesh::getDistanceNumerically(MWMath::Point3D pGlobal, bool signedDistance) {
    // 1. In lokales System transformieren
    MWMath::Point3D diff = pGlobal - PositionGlobal;
    MWMath::Point3D pLoc = OrientationGlobal.transposed().transform(diff);

    // R ist der Major-Radius (Mitte bis Ringzentrum)
    // r ist der Minor-Radius (Dicke des Schlauchs)
    
    // Abstand zur Z-Achse in der XY-Ebene
    double dist_xy = std::sqrt(pLoc.x * pLoc.x + pLoc.y * pLoc.y);
    
    // Abstand zum "Kern-Ring"
    double delta_R = dist_xy - R;
    
    // Euklidischer Abstand zur Ring-Mittellinie
    double dist_to_center_line = std::sqrt(delta_R * delta_R + pLoc.z * pLoc.z);
    
    // Finale Distanz zur Oberfläche
    return dist_to_center_line - r;
}

casadi::MX SSTorusMesh::constraintJacobian(casadi::MX gamma, casadi::MX q) {
    casadi::MX phi = q(casadi::Slice(0, 3));
    casadi::MX d1  = q(casadi::Slice(3, 6));
    casadi::MX d2  = q(casadi::Slice(6, 9));
    casadi::MX d3  = q(casadi::Slice(9, 12));
    casadi::MX diff = gamma - phi;

    casadi::MX dot1 = casadi::MX::dot(diff, d1);
    casadi::MX dot2 = casadi::MX::dot(diff, d2);
    casadi::MX dot3 = casadi::MX::dot(diff, d3);

    casadi::MX term1 = pow(dot1, 2) + pow(dot2, 2) + pow(dot3, 2) + pow(R, 2) - pow(r, 2);
    
    // grad = 2 * term1 * (2*dot1*d1 + 2*dot2*d2 + 2*dot3*d3) - 8*a^2*(dot1*d1 + dot2*d2)
    // (Hinweis: Deine MATLAB Vorlage hat 2*dot1*d1... das entspricht der Ableitung der inneren Klammer)
    return 2.0 * term1 * (2.0 * dot1 * d1 + 2.0 * dot2 * d2 + 2.0 * dot3 * d3) - 
           8.0 * pow(R, 2) * (dot1 * d1 + dot2 * d2);
}

*/
// SCALED TORUS
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
    casadi::MX phi = q(casadi::Slice(0, 3));
    casadi::MX d1  = q(casadi::Slice(3, 6));
    casadi::MX d2  = q(casadi::Slice(6, 9));
    casadi::MX d3  = q(casadi::Slice(9, 12));
    casadi::MX diff = gamma - phi;

    // Skalierte Projektionen
    casadi::MX dot1 = casadi::MX::dot(diff, d1) / A;
    casadi::MX dot2 = casadi::MX::dot(diff, d2) / B;
    casadi::MX dot3 = casadi::MX::dot(diff, d3) / C;

    casadi::MX term1 = pow(dot1, 2) + pow(dot2, 2) + pow(dot3, 2) + pow(R, 2) - pow(r, 2);
    
    // Ableitungen der skalierten Komponenten: d/d_gamma (dot1) = d1 / A
    casadi::MX der_dot1 = (2.0 * dot1 / A) * d1;
    casadi::MX der_dot2 = (2.0 * dot2 / B) * d2;
    casadi::MX der_dot3 = (2.0 * dot3 / C) * d3;

    return 2.0 * term1 * (der_dot1 + der_dot2 + der_dot3) - 
           8.0 * pow(R, 2) * (der_dot1 + der_dot2);
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

////// UNSCALED TORUS //////
/* 
casadi::MX SSEllipsoidMesh::constraintDistance(casadi::MX gamma, casadi::MX q) {
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

casadi::MX SSEllipsoidMesh::constraintJacobian(casadi::MX gamma, casadi::MX q) {
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
 */

void SSMesh::InitializeMesh()
{
    if (Parent){
        Parent->Meshes.push_back(shared_from_this());
        qDebug() << "Mesh " << QString::fromStdString(Name) << " added to Parent " << QString::fromStdString(Parent->Name);
        
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
        OrientationGlobal = Parent->OrientationGlobal * Orientation2ParentRel;
        PositionGlobal = Parent->PositionGlobal + OrientationGlobal.transform(Position2ParentRelInParentFrame);
        if (Name == "_"){
            qDebug() << "    " << this->Name.c_str() << ": Pos: " << QString::fromStdString(PositionGlobal.print()) << ", Ori: " << QString::fromStdString(OrientationGlobal.print());
        }
    } else {
        PositionGlobal = Position2ParentRelInParentFrame; // Falls kein Parent, dann ist die lokale Position die globale
        OrientationGlobal = Orientation2ParentRel; // Falls kein Parent, dann ist die lokale Rotation die globale
    }

    
};
