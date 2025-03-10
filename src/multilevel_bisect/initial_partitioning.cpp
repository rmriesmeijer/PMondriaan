#include <limits>
#include <string>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/KLFM/KLFM.hpp"
#include "multilevel_bisect/initial_partitioning.hpp"
#include "multilevel_bisect/label_propagation.hpp"

namespace pmondriaan {

/**
 * Creates an initial partitioning for hypergraph H. Returns the cutsize of the solution found.
 */
long initial_partitioning(pmondriaan::hypergraph& H,
                          long max_weight_0,
                          long max_weight_1,
                          pmondriaan::options& opts,
                          std::mt19937& rng,
                          std::string breaking_mode) {
    // Breaking up edges of size 3 if required and the majority of hyperedges are size 2 or 3.
    int e23count = -1;
    if (breaking_mode == "break_triples_in_initial_partitioning") {
        //pmondriaan::interval labels = {0,1};
        //bisect_random(H, max_weight_0, max_weight_1, 0, H.size(), labels, rng);
        e23count = 0;
        for (auto n = 0u; n < H.nets().size(); n++) {
            if(H.nets()[n].size() < 4) {
                e23count += 1;
            }
            else {
                e23count -= 2;
            }
        }

        if(e23count < 0)
            simplify_duplicate_nets(H);
        else
            break_triples(H);
    }

    if (breaking_mode == "none") {
        simplify_duplicate_nets(H);
    }

    auto L_best = std::vector<long>(H.size());
    long best_cut = std::numeric_limits<long>::max();
    long best_imbalance = std::numeric_limits<long>::max();
    auto time = bulk::util::timer();
    //std::cout << "Nets size " << H.nets().size();
    for (long i = 0; i < 10; i++) {
        time.get();
        // counts of all labels for each net
        auto C =
        std::vector<std::vector<long>>(H.nets().size(), std::vector<long>(2, 0));

        auto L = label_propagation_bisect(H, C, opts.lp_max_iterations,
                                          max_weight_0, max_weight_1, rng);
        for (auto i = 0u; i < H.size(); i++) {
            H(i).set_part(L[i]);
        }

        // std::cout << "time lp: " << time.get_change() << "(round " << i << ")\n";

        auto cut = pmondriaan::KLFM(H, C, H.weight_part(0), H.weight_part(1),
                                    max_weight_0, max_weight_1, opts, rng);

        cut = pmondriaan::KLFM(H, C, H.weight_part(0), H.weight_part(1),
                                    max_weight_0, max_weight_1, opts, rng);

        // std::cout << "time KLFM: " << time.get_change() << "(round " << i << ")\n";

        long imbalance =
        std::max(H.weight_part(0) - max_weight_0, H.weight_part(1) - max_weight_1);

        if (((cut < best_cut) && (imbalance <= 0)) ||
            (((cut == best_cut) || (best_imbalance > 0)) && (imbalance < best_imbalance))) {
            for (auto i = 0u; i < H.size(); i++) {
                L_best[i] = H(i).part();
            }
            best_cut = cut;
            best_imbalance = imbalance;
        }
    }

    for (auto i = 0u; i < H.size(); i++) {
        H(i).set_part(L_best[i]);
    }

    if (breaking_mode == "break_triples_in_initial_partitioning" && e23count >= 0) {
        return best_cut / 2;
    }

    return best_cut;
}


} // namespace pmondriaan
