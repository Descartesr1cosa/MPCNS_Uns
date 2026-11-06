// MarsUns <setup-file...> [name=value ...]     ten-ion Mars Hall-MHD (run with mpirun, one rank per part file)
#include "app/App.h"
#include "driver/MarsSolver.h"

int main(int argc, char **argv)
{
    return mhd::run_app(argc, argv, "MarsUns", [](MPI_Comm comm, const uns::LocalMesh &lm, const uns::Params &p) {
        mars::MarsSolver solver(comm, lm, mars::MarsParams::from_params(p));
        solver.initialize();
        solver.run();
    });
}
