#pragma once
// Named fields stored per entity number: value(i, comp) = data[i*ncomp + comp].
//  node  : geometry / scalar potentials
//  edge  : 1-cochains  (E, E_hall, J)         [line integrals]
//  face  : 2-cochains  (B)                    [fluxes]
//  cell  : cell averages (U_H, U_Na, reconstructed B, J)
#include <map>
#include <string>
#include <vector>

#include "uns/complex/CellComplex.h"

namespace mercury
{
using uns::Dim;

struct Field
{
    Dim dim;
    int ncomp;
    std::vector<double> data;
    double *ptr() { return data.data(); }
};

class FieldSet
{
public:
    explicit FieldSet(const uns::CellComplex &cx) : cx_(cx) {}
    Field &add(const std::string &name, Dim dim, int ncomp = 1);
    Field &get(const std::string &name);
    bool has(const std::string &name) const { return f_.count(name) != 0; }
    const std::map<std::string, Field> &all() const { return f_; }
private:
    const uns::CellComplex &cx_;
    std::map<std::string, Field> f_;
};
} // namespace mercury
