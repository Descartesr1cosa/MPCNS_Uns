#pragma once
// MPI execution of a CommPlan on entity-numbered arrays (layout: id*ncomp + comp).
//
//   forward(d)        owner value  -> alias copies   (alias = alias_sign * owner)
//   reduce_edges()    alias partial sums -> owner    (owner += alias_sign * partial)
//   assemble_edges()  reduce_edges + forward: every copy of an edge ends with the same
//                     total.  This is the one assembly path shared by ideal, Hall,
//                     resistive and artificial EMF contributions.
#include <mpi.h>

#include <vector>

#include "uns/parallel/LocalMesh.h"

namespace uns::comm
{
class Exchanger
{
public:
    Exchanger(MPI_Comm comm, const LocalMesh &lm);
    void forward(Dim d, double *data, int ncomp = 1);
    void reduce_edges(double *data, int ncomp = 1);
    void assemble_edges(double *data, int ncomp = 1);
    // timed collectives: the time spent here and in the exchanges above is "communication" (including waiting)
    void allreduce(double *v, int n, MPI_Op op);
    double comm_seconds() const { return comm_s_; }
    void reset_comm_clock() { comm_s_ = 0; }
private:
    double comm_s_ = 0;
    void exchange_(const CommLists &send, const CommLists &recv, const EntityInfo &ent, double *data, int ncomp, bool accumulate,
                   bool sign_on_send);
    MPI_Comm comm_;
    const LocalMesh &lm_;
    std::vector<double> sbuf_, rbuf_;
};
} // namespace uns::comm
