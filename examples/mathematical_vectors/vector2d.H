#ifndef VECTOR_2D_H
#define VECTOR_2D_H

#include <iostream>
#include <cmath>

class Vector2d {

private:

    // our member data

    double x;
    double y;

public:

    // default constructor

    Vector2d()
        : x{0.0}, y{0.0}
    {}

    // another constructor

    Vector2d(double _x, double _y)
        : x{_x}, y{_y}
    {}

    // setters

    inline void set_x(double _x) {x = _x;}

    inline void set_y(double _y) {y = _y;}

    // operators

    // add two vectors

    Vector2d operator+ (const Vector2d& vec) {
        Vector2d vec_out(x + vec.x, y + vec.y);
        return vec_out;
    }

    // subtract two vectors

    Vector2d operator- (const Vector2d& vec) {
        Vector2d vec_out(x - vec.x, y - vec.y);
        return vec_out;
    }

    // unary minus

    Vector2d operator- () {
        Vector2d vec_out(-x, -y);
        return vec_out;
    }

    // multiplication

    // vec * a

    Vector2d operator* (double a) {
        Vector2d vec_out(a * x, a * y);
        return vec_out;
    }

    // a * vec -- needs to be a friend

    friend Vector2d operator* (double a, const Vector2d& vec);

    // division by a scalar

    Vector2d operator/ (double a) {
        Vector2d vec_out(x / a, y / a);
        return vec_out;
    }

    // compound operators

    Vector2d& operator+= (const Vector2d& vec) {
        x += vec.x;
        y += vec.y;
        return *this;
    }

    Vector2d& operator-= (const Vector2d& vec) {
        x -= vec.x;
        y -= vec.y;
        return *this;
    }

    Vector2d& operator*= (double a) {
        x *= a;
        y *= a;
        return *this;
    }

    Vector2d& operator/= (double a) {
        x /= a;
        y /= a;
        return *this;
    }

    // comparison operators

    bool operator== (const Vector2d& vec) const {
        return x == vec.x && y == vec.y;
    }

    bool operator!= (const Vector2d& vec) const {
        return x != vec.x || y != vec.y;
    }

    double abs() {
        return std::sqrt(x * x + y * y);
    }

    double dot(const Vector2d& vec) {
        return x * vec.x + y * vec.y;
    }

    // << and >> are not class members, but need access to the member data

    friend std::ostream& operator<< (std::ostream& os, const Vector2d& v);
    friend std::istream& operator>> (std::istream& is, Vector2d& v);
};

inline
std::ostream& operator<< (std::ostream& os, const Vector2d& vec)
{
    os << "(" << vec.x << ", " << vec.y << ")";
    return os;
}

inline
std::istream& operator>> (std::istream& is, Vector2d& vec)
{
    is >> vec.x >> vec.y;
    return is;
}

inline
Vector2d operator*(double a, const Vector2d& vec) {
        Vector2d vec_out(a * vec.x, a * vec.y);
        return vec_out;
}

#endif
