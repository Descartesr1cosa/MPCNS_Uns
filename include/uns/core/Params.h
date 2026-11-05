#pragma once
// Key/value setup file in the legacy style:  "<type> <name> <value>   // comment"
// (type is one of int/double/string/bool and is only a hint).  Command-line
// overrides use  name=value.
#include <map>
#include <string>
#include <vector>

namespace uns
{
class Params
{
public:
    void load_file(const std::string &path);
    void set_override(const std::string &assignment); // "name=value"
    bool has(const std::string &k) const { return kv_.count(k) != 0; }
    int get_int(const std::string &k, int def) const;
    double get_double(const std::string &k, double def) const;
    bool get_bool(const std::string &k, bool def) const;
    std::string get_string(const std::string &k, const std::string &def) const;
    std::vector<double> get_doubles(const std::string &k) const; // whitespace/comma separated
    std::vector<std::string> keys() const;
private:
    std::map<std::string, std::string> kv_;
};
} // namespace uns
