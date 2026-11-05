#include "hall/HallEMF.h"

#include <limits>

namespace mercury
{
void add_hall_emf(const HallConfig &cfg, EdgeEMF &, double *)
{
    (void)cfg; // TODO(Hall): see docs/HallTodo.md
}
double hall_time_step_limit(const HallConfig &) { return std::numeric_limits<double>::infinity(); }
} // namespace mercury
