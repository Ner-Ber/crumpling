#pragma once
#include <cmath>
#include <algorithm>

namespace v3_math {

inline double msqrt(double x) { return std::sqrt(x); }

template <typename T>
inline T max(T a, T b) {
    return (std::max)(a, b);
}

template <typename T>
inline T min(T a, T b) {
    return (std::min)(a, b);
}

}  // namespace v3_math
