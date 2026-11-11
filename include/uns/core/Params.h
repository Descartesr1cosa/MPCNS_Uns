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
    void load_file(const std::string &path);          // legacy "type name value" file
    void load_yaml(const std::string &path);          // flat "key: value" file (paths / launch settings); see load_args
    // Command line = any mix of config files (*.yml / *.yaml / other = legacy txt, later wins) and name=value overrides.
    // A yml file may list further setting files with  settings: [a.txt, b.txt]  (loaded first, the yml keys win).
    // With chdir_to_yml the process enters the directory of the first yml file, so every relative path is relative to the case folder.
    // With no file argument (nothing, or only name=value), ./config.yml is used when it exists.  Returns false when nothing was loaded.
    bool load_args(const std::vector<std::string> &args, bool chdir_to_yml = true);
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
