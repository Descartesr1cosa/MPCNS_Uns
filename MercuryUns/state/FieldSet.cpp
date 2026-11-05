#include "state/FieldSet.h"

#include <stdexcept>

namespace mercury
{
Field &FieldSet::add(const std::string &name, Dim dim, int ncomp)
{
    if (f_.count(name)) throw std::runtime_error("field '" + name + "' already exists");
    Field f{dim, ncomp, std::vector<double>(static_cast<std::size_t>(cx_.count(dim)) * ncomp, 0.0)};
    return f_.emplace(name, std::move(f)).first->second;
}
Field &FieldSet::get(const std::string &name)
{
    auto it = f_.find(name);
    if (it == f_.end()) throw std::runtime_error("no field '" + name + "'");
    return it->second;
}
} // namespace mercury
