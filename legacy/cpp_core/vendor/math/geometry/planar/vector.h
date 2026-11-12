#pragma once
#include <cmath>

namespace v3_math {
namespace spatial {

template <typename T>
struct vector {
    T x = T();
    T y = T();
    T z = T();

    vector() = default;
    vector(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

    vector& operator+=(const vector& o) {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }
    vector& operator-=(const vector& o) {
        x -= o.x;
        y -= o.y;
        z -= o.z;
        return *this;
    }
    vector& operator*=(T s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    T length2() const { return x * x + y * y + z * z; }
    T length() const { return std::sqrt(length2()); }
};

template <typename T>
inline vector<T> operator+(vector<T> a, const vector<T>& b) {
    return a += b;
}
template <typename T>
inline vector<T> operator-(vector<T> a, const vector<T>& b) {
    return a -= b;
}
template <typename T>
inline vector<T> operator*(vector<T> a, T s) {
    return a *= s;
}
template <typename T>
inline vector<T> operator*(T s, vector<T> a) {
    return a *= s;
}

}  // namespace spatial
}  // namespace v3_math
