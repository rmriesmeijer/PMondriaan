#include <random>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include "util/interval.hpp"

namespace pmondriaan {

/**
 * Creates an initial partitioning for hypergraph H.
 */
void initial_partitioning(bulk::world& world,
                          pmondriaan::hypergraph& H,
                          long max_weight_0,
                          long max_weight_1,
                          interval labels,
                          std::mt19937& rng);

} // namespace pmondriaan
