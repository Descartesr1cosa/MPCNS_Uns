#pragma once
#include <string>
#include <vector>

#include "uns/partition/CellGraph.h"

namespace uns
{
enum class PartMethod
{
    Morton,  // space-filling curve, weighted prefix cut
    Hilbert, // space-filling curve, weighted prefix cut
    Rcb,     // weighted recursive coordinate bisection
    Metis    // METIS k-way (only when built with UNS_ENABLE_METIS)
};
PartMethod parse_method(const std::string &s);
std::string to_string(PartMethod m);

struct PartOptions
{
    PartMethod method = PartMethod::Hilbert;
    int refine_passes = 8;     // greedy boundary refinement (0 = off)
    double tolerance = 0.03;   // allowed load imbalance during refinement
};

// Returns the rank of every cell.
std::vector<int> partition_cells(const CellGraph &g, int nparts, const PartOptions &opt);

struct PartitionReport
{
    int nparts = 0;
    std::vector<double> load;
    std::vector<Id> cells;
    std::vector<Id> cut_faces;     // faces shared with other ranks
    std::vector<double> cut_weight;// communication weight to other ranks
    std::vector<int> neighbours;
    double imbalance = 0;          // max load / mean load
    double total_cut_weight = 0;
    long total_cut_faces = 0;
    std::string format(const std::string &title) const;
};
PartitionReport evaluate_partition(const CellGraph &g, const std::vector<int> &part, int nparts);
} // namespace uns
