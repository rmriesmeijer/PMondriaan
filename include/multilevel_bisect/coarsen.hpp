#include <string>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "bisect.hpp"
#include "hypergraph/contraction.hpp"
#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Coarses the hypergraph H and returns a hypergraph HC.
 */
pmondriaan::hypergraph coarsen_hypergraph(bulk::world& world,
                                          pmondriaan::hypergraph& H,
                                          pmondriaan::contraction& C,
                                          pmondriaan::options& opts,
                                          std::string sampling_mode);


/**
 * Sends match request to the owners of the best matches found using the
 * improduct computation. Returns the local matches.
 */
void request_matches(pmondriaan::hypergraph& H,
                     pmondriaan::contraction& C,
                     bulk::queue<int, long, int[]>& sample_queue,
                     bulk::queue<int, int>& accepted_matches,
                     const std::vector<int>& indices_samples,
                     pmondriaan::options& opts);

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           const std::vector<int> samples,
                                           bulk::queue<int, long, int[]>& matches,
                                           std::vector<bool>& matched);


} // namespace pmondriaan
