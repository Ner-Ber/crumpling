#pragma once
#include <cstdlib>
#include <string>

namespace v3_general {

inline int stoi(const std::string& s) {
    return static_cast<int>(std::strtol(s.c_str(), nullptr, 10));
}

inline double stod(const std::string& s) {
    return std::strtod(s.c_str(), nullptr);
}

}  // namespace v3_general
