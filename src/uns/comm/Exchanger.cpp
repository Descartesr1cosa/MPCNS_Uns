#include "uns/comm/Exchanger.h"

namespace uns::comm
{
Exchanger::Exchanger(MPI_Comm comm, const LocalMesh &lm) : comm_(comm), lm_(lm) {}

// send list refers to entities whose values leave this rank, recv list to the ones that arrive.
// forward : send = owned (value as is), recv = ghost (value = alias_sign * received)
// reduce  : send = alias partials (value * alias_sign), recv = owned (accumulate)
void Exchanger::exchange_(const CommLists &send, const CommLists &recv, const EntityInfo &ent, double *data, int nc,
                          bool accumulate, bool sign_on_send)
{
    sbuf_.resize(send.total() * nc);
    rbuf_.resize(recv.total() * nc);
    std::vector<MPI_Request> req;
    req.reserve(send.peers.size() + recv.peers.size());
    for (std::size_t i = 0; i < recv.peers.size(); ++i)
    {
        req.emplace_back();
        MPI_Irecv(rbuf_.data() + recv.ptr[i] * nc, static_cast<int>((recv.ptr[i + 1] - recv.ptr[i]) * nc), MPI_DOUBLE,
                  recv.peers[i], 77, comm_, &req.back());
    }
    for (std::size_t i = 0; i < send.peers.size(); ++i)
    {
        for (std::int64_t k = send.ptr[i]; k < send.ptr[i + 1]; ++k)
        {
            const Id id = send.idx[k];
            const double s = sign_on_send ? ent.alias_sign[id] : 1.0;
            for (int c = 0; c < nc; ++c) sbuf_[k * nc + c] = s * data[static_cast<std::size_t>(id) * nc + c];
        }
        req.emplace_back();
        MPI_Isend(sbuf_.data() + send.ptr[i] * nc, static_cast<int>((send.ptr[i + 1] - send.ptr[i]) * nc), MPI_DOUBLE,
                  send.peers[i], 77, comm_, &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    for (std::size_t k = 0; k < recv.idx.size(); ++k)
    {
        const Id id = recv.idx[k];
        const double s = sign_on_send ? 1.0 : ent.alias_sign[id];
        for (int c = 0; c < nc; ++c)
        {
            double &dst = data[static_cast<std::size_t>(id) * nc + c];
            dst = accumulate ? dst + s * rbuf_[k * nc + c] : s * rbuf_[k * nc + c];
        }
    }
}

void Exchanger::forward(Dim d, double *data, int nc)
{
    exchange_(lm_.plan.send[di(d)], lm_.plan.recv[di(d)], lm_.mesh.cx.ent[di(d)], data, nc, false, false);
}

void Exchanger::reduce_edges(double *data, int nc)
{
    exchange_(lm_.plan.reduce_send, lm_.plan.reduce_recv, lm_.mesh.cx.ent[1], data, nc, true, true);
}

void Exchanger::assemble_edges(double *data, int nc)
{
    reduce_edges(data, nc);
    forward(Dim::Edge, data, nc);
}
} // namespace uns::comm
