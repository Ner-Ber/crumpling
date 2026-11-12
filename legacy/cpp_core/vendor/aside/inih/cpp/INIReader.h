#pragma once
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

// Minimal INIReader compatible with the subset used by lmp_swapstat4_main.
class INIReader {
public:
    explicit INIReader(const std::string& filename) {
        std::ifstream in(filename.c_str());
        if (!in) {
            error_ = -1;
            return;
        }
        std::string section;
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0] == ';' || line[0] == '#') continue;
            if (line[0] == '[') {
                auto end = line.find(']');
                if (end == std::string::npos) {
                    error_ = 1;
                    continue;
                }
                section = to_lower(trim(line.substr(1, end - 1)));
                continue;
            }
            auto eq = line.find('=');
            if (eq == std::string::npos) eq = line.find(':');
            if (eq == std::string::npos) {
                error_ = 1;
                continue;
            }
            std::string name = to_lower(trim(line.substr(0, eq)));
            std::string value = trim(line.substr(eq + 1));
            values_[{section, name}] = value;
        }
    }

    int ParseError() const { return error_; }

    double GetReal(const std::string& section, const std::string& name, double def) const {
        auto it = values_.find({to_lower(section), to_lower(name)});
        if (it == values_.end()) return def;
        return std::strtod(it->second.c_str(), nullptr);
    }

    long GetInteger(const std::string& section, const std::string& name, long def) const {
        auto it = values_.find({to_lower(section), to_lower(name)});
        if (it == values_.end()) return def;
        return std::strtol(it->second.c_str(), nullptr, 0);
    }

    bool GetBoolean(const std::string& section, const std::string& name, bool def) const {
        auto it = values_.find({to_lower(section), to_lower(name)});
        if (it == values_.end()) return def;
        std::string v = to_lower(it->second);
        return v == "true" || v == "yes" || v == "on" || v == "1";
    }

private:
    int error_ = 0;
    std::map<std::pair<std::string, std::string>, std::string> values_;

    static std::string trim(std::string s) {
        while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
        while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
        return s;
    }
    static std::string to_lower(std::string s) {
        for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    }
};
