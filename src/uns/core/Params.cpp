#include "uns/core/Params.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace uns
{
namespace
{
std::string strip(std::string s)
{
    const auto c = s.find("//");
    if (c != std::string::npos) s.erase(c);
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    return s.substr(b, s.find_last_not_of(" \t\r\n") - b + 1);
}
bool is_type(const std::string &t) { return t == "int" || t == "double" || t == "string" || t == "bool"; }
} // namespace

void Params::load_file(const std::string &path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open setup file '" + path + "'");
    bool block = false;
    for (std::string line; std::getline(in, line);)
    {
        if (line.find("/*") != std::string::npos) block = true;
        if (block) { if (line.find("*/") != std::string::npos) block = false; continue; }
        line = strip(line);
        if (line.empty()) continue;
        std::istringstream ls(line);
        std::string a, b;
        ls >> a;
        if (is_type(a)) ls >> b; else b = a;
        std::string rest;
        std::getline(ls, rest);
        kv_[b] = strip(rest);
    }
}

void Params::set_override(const std::string &a)
{
    const auto eq = a.find('=');
    if (eq == std::string::npos) throw std::runtime_error("override must be name=value: '" + a + "'");
    kv_[strip(a.substr(0, eq))] = strip(a.substr(eq + 1));
}

std::string Params::get_string(const std::string &k, const std::string &d) const
{
    auto it = kv_.find(k);
    if (it == kv_.end()) return d;
    const std::string &v = it->second;
    if (v.size() >= 2 && v.front() == '"' && v.back() == '"') return v.substr(1, v.size() - 2); // "" = empty string
    return v;
}
int Params::get_int(const std::string &k, int d) const
{
    auto it = kv_.find(k);
    return it == kv_.end() ? d : std::stoi(it->second);
}
double Params::get_double(const std::string &k, double d) const
{
    auto it = kv_.find(k);
    return it == kv_.end() ? d : std::stod(it->second);
}
bool Params::get_bool(const std::string &k, bool d) const
{
    auto it = kv_.find(k);
    if (it == kv_.end()) return d;
    return it->second == "1" || it->second == "true" || it->second == "yes";
}
std::vector<double> Params::get_doubles(const std::string &k) const
{
    std::vector<double> v;
    auto it = kv_.find(k);
    if (it == kv_.end()) return v;
    std::string s = it->second;
    for (char &c : s) if (c == ',') c = ' ';
    std::istringstream is(s);
    for (double x; is >> x;) v.push_back(x);
    return v;
}
std::vector<std::string> Params::keys() const
{
    std::vector<std::string> v;
    for (const auto &kv : kv_) v.push_back(kv.first);
    return v;
}
} // namespace uns
