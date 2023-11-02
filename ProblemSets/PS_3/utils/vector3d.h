#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>

class Vector3D {
public:
    // Default constructor
    double x, y, z;

    Vector3D() : x(0.0), y(0.0), z(0.0) {}

    // Parameterized constructor
    Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Compute the magnitude (norm) of the vector
    double norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }
};

#endif // VECTOR3D_H
