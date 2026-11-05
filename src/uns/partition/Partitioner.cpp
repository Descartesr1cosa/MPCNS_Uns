#include "uns/partition/Partitioner.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace uns
{
PartMethod parse_method(const std::string &s)
{
    if (s == "morton") return PartMethod::Morton;
    if (s == "hilbert" || s == "sfc") return PartMethod::Hilbert;
    if (s == "rcb") return PartMethod::Rcb;
    if (s == "metis") return PartMethod::Metis;
    throw std::runtime_error("unknown partition method '" + s + "' (morton|hilbert|rcb|metis)");
}
std::string to_string(PartMethod m)
{
    switch (m)
    {
    case PartMethod::Morton: return "morton";
    case PartMethod::Hilbert: return "hilbert";
    case PartMethod::Rcb: return "rcb";
    case PartMethod::Metis: return "metis";
    }
    return "?";
}

namespace
{
constexpr int kBits = 20;

std::uint64_t interleave3(std::uint32_t x, std::uint32_t y, std::uint32_t z)
{
    std::uint64_t k = 0;
    for (int b = kBits - 1; b >= 0; --b)
        k = (k << 3) | (((x >> b) & 1u) << 2) | (((y >> b) & 1u) << 1) | ((z >> b) & 1u);
    return k;
}

// Skilling's transpose algorithm: axes -> Hilbert index.
std::uint64_t hilbert3(std::uint32_t X[3])
{
    const std::uint32_t M = 1u << (kBits - 1);
    for (std::uint32_t Q = M; Q > 1; Q >>= 1)
    {
        const std::uint32_t P = Q - 1;
        for (int i = 0; i < 3; ++i)
        {
            if (X[i] & Q) X[0] ^= P;
            else { const std::uint32_t t = (X[0] ^ X[i]) & P; X[0] ^= t; X[i] ^= t; }
        }
    }
    for (int i = 1; i < 3; ++i) X[i] ^= X[i - 1];
    std::uint32_t t = 0;
    for (std::uint32_t Q = M; Q > 1; Q >>= 1)
        if (X[2] & Q) t ^= Q - 1;
    for (int i = 0; i < 3; ++i) X[i] ^= t;
    return interleave3(X[0], X[1], X[2]);
}

void bbox(const CellGraph &g, Vec3 &lo, Vec3 &hi)
{
    lo = hi = g.centre.empty() ? Vec3() : g.centre[0];
    for (const Vec3 &p : g.centre)
    {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
}

std::vector<int> cut_by_weight(const CellGraph &g, const std::vector<Id> &order, int nparts)
{
    std::vector<int> part(g.n(), 0);
    const double total = std::accumulate(g.vw.begin(), g.vw.end(), 0.0);
    double acc = 0;
    for (Id c : order)
    {
        const double mid = acc + 0.5 * g.vw[c]; // assign by the cell's weight midpoint
        part[c] = std::min(nparts - 1, static_cast<int>(mid / total * nparts));
        acc += g.vw[c];
    }
    return part;
}

std::vector<int> sfc(const CellGraph &g, int nparts, bool hilbert)
{
    Vec3 lo, hi;
    bbox(g, lo, hi);
    const double ext = std::max({hi.x - lo.x, hi.y - lo.y, hi.z - lo.z, 1e-300});
    const double scale = ((1u << kBits) - 1) / ext;
    std::vector<std::uint64_t> key(g.n());
    for (Id c = 0; c < g.n(); ++c)
    {
        std::uint32_t q[3] = {static_cast<std::uint32_t>((g.centre[c].x - lo.x) * scale),
                              static_cast<std::uint32_t>((g.centre[c].y - lo.y) * scale),
                              static_cast<std::uint32_t>((g.centre[c].z - lo.z) * scale)};
        key[c] = hilbert ? hilbert3(q) : interleave3(q[0], q[1], q[2]);
    }
    std::vector<Id> order(g.n());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](Id a, Id b) { return key[a] != key[b] ? key[a] < key[b] : a < b; });
    return cut_by_weight(g, order, nparts);
}

void rcb_rec(const CellGraph &g, std::vector<Id> &ids, std::size_t b, std::size_t e, int p0, int np, std::vector<int> &part)
{
    if (np == 1)
    {
        for (std::size_t i = b; i < e; ++i) part[ids[i]] = p0;
        return;
    }
    Vec3 lo = g.centre[ids[b]], hi = lo;
    for (std::size_t i = b; i < e; ++i)
    {
        const Vec3 &p = g.centre[ids[i]];
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    const double dx = hi.x - lo.x, dy = hi.y - lo.y, dz = hi.z - lo.z;
    const int axis = (dx >= dy && dx >= dz) ? 0 : (dy >= dz ? 1 : 2);
    auto coord = [&](Id c) { return axis == 0 ? g.centre[c].x : axis == 1 ? g.centre[c].y : g.centre[c].z; };
    std::sort(ids.begin() + b, ids.begin() + e, [&](Id a, Id c) { return coord(a) != coord(c) ? coord(a) < coord(c) : a < c; });
    const int npl = np / 2;
    double total = 0;
    for (std::size_t i = b; i < e; ++i) total += g.vw[ids[i]];
    const double target = total * npl / np;
    double acc = 0;
    std::size_t m = b;
    while (m < e - 1 && acc + 0.5 * g.vw[ids[m]] < target) acc += g.vw[ids[m++]];
    if (m == b) m = b + 1;
    rcb_rec(g, ids, b, m, p0, npl, part);
    rcb_rec(g, ids, m, e, p0 + npl, np - npl, part);
}

// Greedy boundary refinement: move boundary cells to the neighbouring part with
// the largest cut gain while keeping every load within (1+tol) of the mean.
void refine(const CellGraph &g, std::vector<int> &part, int nparts, int passes, double tol)
{
    std::vector<double> load(nparts, 0.0);
    for (Id c = 0; c < g.n(); ++c) load[part[c]] += g.vw[c];
    const double mean = std::accumulate(load.begin(), load.end(), 0.0) / nparts;
    const double cap = mean * (1.0 + tol), floor_ = mean * (1.0 - tol);
    std::vector<double> conn(nparts);
    for (int pass = 0; pass < passes; ++pass)
    {
        long moved = 0;
        for (Id c = 0; c < g.n(); ++c)
        {
            const int p = part[c];
            std::fill(conn.begin(), conn.end(), 0.0);
            bool boundary = false;
            for (std::int64_t k = g.xadj[c]; k < g.xadj[c + 1]; ++k)
            {
                conn[part[g.adj[k]]] += g.adjw[k];
                boundary |= part[g.adj[k]] != p;
            }
            if (!boundary) continue;
            int best = p;
            double bestgain = 1e-12;
            for (int q = 0; q < nparts; ++q)
            {
                if (q == p || conn[q] <= 0) continue;
                const double gain = conn[q] - conn[p];
                if (gain > bestgain && load[q] + g.vw[c] <= cap && load[p] - g.vw[c] >= floor_)
                {
                    best = q;
                    bestgain = gain;
                }
            }
            if (best != p)
            {
                load[p] -= g.vw[c];
                load[best] += g.vw[c];
                part[c] = best;
                ++moved;
            }
        }
        if (moved == 0) break;
    }
}
} // namespace

std::vector<int> partition_cells(const CellGraph &g, int nparts, const PartOptions &opt)
{
    if (nparts < 1) throw std::runtime_error("partition_cells: nparts < 1");
    if (nparts > g.n()) throw std::runtime_error("partition_cells: more parts than cells");
    std::vector<int> part;
    switch (opt.method)
    {
    case PartMethod::Morton: part = sfc(g, nparts, false); break;
    case PartMethod::Hilbert: part = sfc(g, nparts, true); break;
    case PartMethod::Rcb:
    {
        std::vector<Id> ids(g.n());
        std::iota(ids.begin(), ids.end(), 0);
        part.assign(g.n(), 0);
        rcb_rec(g, ids, 0, ids.size(), 0, nparts, part);
        break;
    }
    case PartMethod::Metis:
        throw std::runtime_error("partition method 'metis' requested but this build has no METIS "
                                 "(configure with -DUNS_ENABLE_METIS=ON and a METIS install)");
    }
    if (opt.refine_passes > 0 && nparts > 1) refine(g, part, nparts, opt.refine_passes, opt.tolerance);
    return part;
}

PartitionReport evaluate_partition(const CellGraph &g, const std::vector<int> &part, int nparts)
{
    PartitionReport r;
    r.nparts = nparts;
    r.load.assign(nparts, 0.0);
    r.cells.assign(nparts, 0);
    r.cut_faces.assign(nparts, 0);
    r.cut_weight.assign(nparts, 0.0);
    r.neighbours.assign(nparts, 0);
    std::vector<std::set<int>> nb(nparts);
    for (Id c = 0; c < g.n(); ++c)
    {
        r.load[part[c]] += g.vw[c];
        ++r.cells[part[c]];
        for (std::int64_t k = g.xadj[c]; k < g.xadj[c + 1]; ++k)
            if (part[g.adj[k]] != part[c])
            {
                ++r.cut_faces[part[c]];
                r.cut_weight[part[c]] += g.adjw[k];
                nb[part[c]].insert(part[g.adj[k]]);
            }
    }
    double sum = 0, mx = 0;
    for (int p = 0; p < nparts; ++p)
    {
        sum += r.load[p];
        mx = std::max(mx, r.load[p]);
        r.neighbours[p] = static_cast<int>(nb[p].size());
        r.total_cut_faces += r.cut_faces[p];
        r.total_cut_weight += r.cut_weight[p];
    }
    r.total_cut_faces /= 2;     // every cut face is counted from both sides
    r.total_cut_weight *= 0.5;
    r.imbalance = mx / (sum / nparts);
    return r;
}

std::string PartitionReport::format(const std::string &title) const
{
    std::ostringstream s;
    s << "=== partition report: " << title << " (" << nparts << " ranks) ===\n";
    s << "rank      cells        load  load/mean  cut_faces  cut_weight  neighbours\n";
    double sum = 0;
    for (double l : load) sum += l;
    const double mean = sum / nparts;
    for (int p = 0; p < nparts; ++p)
    {
        char buf[200];
        std::snprintf(buf, sizeof buf, "%4d  %9d  %10.2f  %9.3f  %9d  %10.2f  %10d\n", p, cells[p], load[p],
                      load[p] / mean, cut_faces[p], cut_weight[p], neighbours[p]);
        s << buf;
    }
    char buf[200];
    std::snprintf(buf, sizeof buf, "load imbalance (max/mean) = %.4f   cut faces = %ld   cut weight = %.2f\n", imbalance,
                  total_cut_faces, total_cut_weight);
    s << buf;
    return s.str();
}
} // namespace uns
