#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include <casadi/casadi.hpp>

namespace MWMath {

    struct Point3D {
        double x;
        double y;
        double z;

        Point3D() : x(0), y(0), z(0) {}
        Point3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

        // Operatoren
        Point3D operator+(const Point3D& other) const {
            return {x + other.x, y + other.y, z + other.z};
        }

        Point3D operator-(const Point3D& other) const {
            return {x - other.x, y - other.y, z - other.z};
        }

        Point3D operator*(double scalar) const {
            return {x * scalar, y * scalar, z * scalar};
        }

        void operator+=(const Point3D& other) {
            x += other.x;
            y += other.y;
            z += other.z;
        }

        Point3D normed() const {
            double len = std::sqrt(x*x + y*y + z*z);
            if (len < 1e-9) return {0, 0, 0};
            return {x / len, y / len, z / len};
        }

        std::string print() const {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(4); // 4 Nachkommastellen
            ss << "(" << x << ", " << y << ", " << z << ")";
            return ss.str();
        }
    };

    // ---------------------------
    // 3x3 Rotationsmatrix
    // ---------------------------
    struct RotMatrix3x3 {
        double m[3][3]; // Zeilen-major

        RotMatrix3x3() {
            // Identität
            m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0;
            m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0;
            m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0;
        }

        // Konstruktor mit 9 Werten
        // R = (x , y, z) -> axis
        RotMatrix3x3(double a00,double a01,double a02,
                    double a10,double a11,double a12,
                    double a20,double a21,double a22)
        {
            m[0][0]=a00; m[0][1]=a01; m[0][2]=a02;
            m[1][0]=a10; m[1][1]=a11; m[1][2]=a12;
            m[2][0]=a20; m[2][1]=a21; m[2][2]=a22;
        }

        RotMatrix3x3 transposed() const {
            return RotMatrix3x3(
                m[0][0], m[1][0], m[2][0],
                m[0][1], m[1][1], m[2][1],
                m[0][2], m[1][2], m[2][2]
            );
        }

        // Matrix * Point3D
        Point3D transform(const Point3D& p) const {
            return {
                m[0][0]*p.x + m[0][1]*p.y + m[0][2]*p.z,
                m[1][0]*p.x + m[1][1]*p.y + m[1][2]*p.z,
                m[2][0]*p.x + m[2][1]*p.y + m[2][2]*p.z
            };
        }

        double determinant() {
            return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        }

        // Matrix * Point3D
        Point3D operator*(const Point3D& p) const {
            return {
                m[0][0]*p.x + m[0][1]*p.y + m[0][2]*p.z,
                m[1][0]*p.x + m[1][1]*p.y + m[1][2]*p.z,
                m[2][0]*p.x + m[2][1]*p.y + m[2][2]*p.z
            };
        }

        // Matrix * Matrix
        RotMatrix3x3 operator*(const RotMatrix3x3& other) const {
            RotMatrix3x3 result;
            for(int i=0;i<3;i++) {
                for(int j=0;j<3;j++) {
                    result.m[i][j] = 0.0;
                    for(int k=0;k<3;k++) {
                        result.m[i][j] += m[i][k] * other.m[k][j];
                    }
                }
            }
            return result;
        }

        RotMatrix3x3& operator*=(const RotMatrix3x3& other) {
        double temp[3][3];
            // 1. Standard Matrix-Multiplikation in temporären Puffer
            // Zeile i, Spalte j
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    temp[i][j] = 0.0;
                    // Skalarprodukt aus Zeile von 'this' und Spalte von 'other'
                    for (int k = 0; k < 3; ++k) {
                        temp[i][j] += m[i][k] * other.m[k][j];
                    }
                }
            }
            // 2. Ergebnis zurück in eigene Matrix kopieren
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    m[i][j] = temp[i][j];
                }
            }
            // Referenz auf sich selbst zurückgeben (für Chaining wie (A *= B) *= C)
            return *this;
        }

        std::string print() const {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(3);
            ss << "\n    [" << m[0][0] << ", " << m[0][1] << ", " << m[0][2] << "]\n"
               << "    [" << m[1][0] << ", " << m[1][1] << ", " << m[1][2] << "]\n"
               << "    [" << m[2][0] << ", " << m[2][1] << ", " << m[2][2] << "]";
            return ss.str();
        }
    };


    inline double distance(const Point3D& a, const Point3D& b) {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        double dz = a.z - b.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }


    inline RotMatrix3x3 axisAngle(const Point3D& axis, double angleDeg) {
        
        double angleRad = angleDeg * (M_PI / 180.0); 

        
        double length = std::sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
        if (length < 1e-9) return RotMatrix3x3(); // Identitätsmatrix bei Null-Achse

        double x = axis.x / length;
        double y = axis.y / length;
        double z = axis.z / length;

        double c = std::cos(angleRad);
        double s = std::sin(angleRad);
        double t = 1.0 - c;

        
        return RotMatrix3x3(
            t*x*x + c,   t*x*y - s*z, t*x*z + s*y,
            t*x*y + s*z, t*y*y + c,   t*y*z - s*x,
            t*x*z - s*y, t*y*z + s*x, t*z*z + c
        );
}

    inline Point3D Rot3x3ToEulerZYX(const RotMatrix3x3& R) {
        Point3D euler;
        if (R.m[2][0] < 1) {
            if (R.m[2][0] > -1) {
                euler.y = std::asin(-R.m[2][0]); // beta
                euler.x = std::atan2(R.m[2][1], R.m[2][2]); // alpha
                euler.z = std::atan2(R.m[1][0], R.m[0][0]); // gamma
            } else { // R.m[2][0] == -1
                euler.y = M_PI / 2;
                euler.x = -std::atan2(-R.m[1][2], R.m[1][1]);
                euler.z = 0;
            }
        } else { // R.m[2][0] == 1
            euler.y = -M_PI / 2;
            euler.x = std::atan2(-R.m[1][2], R.m[1][1]);
            euler.z = 0;
        }
        return euler;
    }

    inline RotMatrix3x3 eulerZYXToRot3x3(double xDeg, double yDeg, double zDeg) {
        // Einfachste Variante: Kombiniere 3 Achsen-Rotationen
        return axisAngle({0,0,1}, zDeg) * axisAngle({0,1,0}, yDeg) * axisAngle({1,0,0}, xDeg);
    }

    inline std::vector<double> AxisAngleFromRotMatrix(const RotMatrix3x3& R) {
        std::vector<double> result(4, 0.0); // x, y, z, angleDeg

        double angleRad = std::acos((R.m[0][0] + R.m[1][1] + R.m[2][2] - 1.0) / 2.0);
        result[3] = angleRad * (180.0 / M_PI); // in Grad

        if (std::abs(angleRad) < 1e-6) {
            result[0] = 1.0; result[1] = 0.0; result[2] = 0.0;
            return result;
        }

        double denom = 2.0 * std::sin(angleRad);
        result[0] = (R.m[2][1] - R.m[1][2]) / denom;
        result[1] = (R.m[0][2] - R.m[2][0]) / denom;
        result[2] = (R.m[1][0] - R.m[0][1]) / denom;

        return result;
    }

    inline std::vector<casadi::MX> MWPoint2Cas(const Point3D& p) {
        std::vector<casadi::MX> casadiPoint(3);
        casadiPoint[0] = casadi::MX(p.x);
        casadiPoint[1] = casadi::MX(p.y);
        casadiPoint[2] = casadi::MX(p.z);
        return casadiPoint;
    }

    inline Point3D Cas2MWPoint(const std::vector<casadi::MX>& casadiPoint) {
        if (casadiPoint.size() != 3) {
            throw std::invalid_argument("Cas2MWPoint: Input vector must have exactly 3 elements.");
        }
        Point3D p;
        p.x = static_cast<double>(casadi::DM(casadiPoint[0]).scalar());
        p.y = static_cast<double>(casadi::DM(casadiPoint[1]).scalar());
        p.z = static_cast<double>(casadi::DM(casadiPoint[2]).scalar());
        return p;
    }

    inline std::vector<casadi::MX> MWPoint2CasConst(const Point3D& p) {
        std::vector<casadi::MX> casadiPoint(3);
        casadiPoint[0] = casadi::MX::vertcat({casadi::MX(p.x)});
        casadiPoint[1] = casadi::MX::vertcat({casadi::MX(p.y)});
        casadiPoint[2] = casadi::MX::vertcat({casadi::MX(p.z)});
        return casadiPoint;
    }

    inline casadi::MX MWPoint2MXConst(const Point3D& p) {
        casadi::MX casadiPoint = casadi::MX::vertcat({casadi::MX(p.x), casadi::MX(p.y), casadi::MX(p.z)});
        return casadiPoint;
    }

    inline std::vector<MWMath::Point3D> MWDoubleVec2Point3DVec(const std::vector<double>& vec) {
        std::vector<MWMath::Point3D> points;
        if (vec.size() % 3 != 0) {
            throw std::invalid_argument("Input vector size must be a multiple of 3.");
        }
        for (size_t i = 0; i < vec.size(); i += 3) {
            points.emplace_back(vec[i], vec[i + 1], vec[i + 2]);
        }
        return points;
    }

    inline casadi::MX MWRotation2MXConst(const RotMatrix3x3& R) {
        casadi::MX R_cas = casadi::MX::zeros(3, 3);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                R_cas(i, j) = R.m[i][j];
            }
        }
        return R_cas;
    }
}


namespace Utility {
    inline void log(const std::string& message, int logLevel = 1) {
        // Einfache Konsolenausgabe
        std::cout << "[LOG]: " << message << std::endl;
    }

    

}