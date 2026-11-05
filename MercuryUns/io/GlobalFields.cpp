#include "io/GlobalFields.h"

#include <algorithm>
#include <numeric>

namespace mercury
{
std::vector<double> gather_entities(MPI_Comm comm, const uns::EntityInfo &ent, const double *values, int nc)
{
    int rank, np;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);
    const int n = ent.n_owned;
    long long mx = -1, gmx = -1;
    for (int i = 0; i < n; ++i) mx = std::max<long long>(mx, ent.gid[i]);
    MPI_Allreduce(&mx, &gmx, 1, MPI_LONG_LONG, MPI_MAX, comm);
    std::vector<int> counts(np), displ(np, 0), vcounts(np), vdispl(np, 0);
    MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
    std::vector<long long> gids(ent.gid.begin(), ent.gid.begin() + n), allg;
    std::vector<double> allv;
    if (rank == 0)
    {
        for (int r = 0; r < np; ++r) { vcounts[r] = counts[r] * nc; if (r) { displ[r] = displ[r - 1] + counts[r - 1]; vdispl[r] = vdispl[r - 1] + vcounts[r - 1]; } }
        allg.resize(std::accumulate(counts.begin(), counts.end(), 0));
        allv.resize(allg.size() * nc);
    }
    MPI_Gatherv(gids.data(), n, MPI_LONG_LONG, allg.data(), counts.data(), displ.data(), MPI_LONG_LONG, 0, comm);
    MPI_Gatherv(values, n * nc, MPI_DOUBLE, allv.data(), vcounts.data(), vdispl.data(), MPI_DOUBLE, 0, comm);
    std::vector<double> out;
    if (rank == 0)
    {
        out.assign(static_cast<std::size_t>(gmx + 1) * nc, 0.0);
        for (std::size_t i = 0; i < allg.size(); ++i)
            for (int k = 0; k < nc; ++k) out[static_cast<std::size_t>(allg[i]) * nc + k] = allv[i * nc + k];
    }
    return out;
}
} // namespace mercury
