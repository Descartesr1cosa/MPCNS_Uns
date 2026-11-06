// MercuryUns <setup-file | -> [name=value ...]     (run with mpirun, one rank per part file)
#include "app/App.h"
#include "driver/Solver.h"

int main(int argc, char **argv)
{
    return mhd::run_app(argc, argv, "MercuryUns", [](MPI_Comm comm, const uns::LocalMesh &lm, const uns::Params &p) {
        mercury::Solver solver(comm, lm, mhd::Config::from_params(p));
        solver.initialize();
        solver.run();
    });
}
