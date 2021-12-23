#include <limits>

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
                          std::mt19937& rng) {

    // pmondriaan::interval labels = {0,1};
    // bisect_random(H, max_weight_0, max_weight_1, 0, H.size(), labels, rng);

    simplify_duplicate_nets(H);

    bool broken = true;

    unsigned int c = 0;
    for (auto n = 0u; n < H.nets().size() && broken == false; n++) {
        if(H.nets()[n].size() == 3) {
            c++;
        }
    }

    auto co = 0;
    for (auto n = 0u; n < H.nets().size(); n++) {
        co += H.nets()[n].cost();
    }
    std::cout << "totco: " << co << "\n";

    if(4 * c > H.nets().size()) {
        break_triples(H);
        broken = true;
    }

    //H.print();
    std::cout << "totco2: " << H.nets().size() << "\n";

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

        auto cs = cutsize(H, opts.metric);
        std::cout << "lp cs: " << cs << "\n";


        // std::cout << "time lp: " << time.get_change() << "(round " << i << ")\n";

        auto cut = pmondriaan::KLFM(H, C, H.weight_part(0), H.weight_part(1),
                                    max_weight_0, max_weight_1, opts, rng);

        std::cout << "kl cs: " << cs << "\n";

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

    if(broken && c != 0) {
        for (auto n = 0u; n < H.nets().size(); n++) {
            long cost = H.nets()[n].cost();
            H.nets()[n].set_cost(cost / 2);
        }
        std::cout << "init: " << best_cut << "\n";

        return best_cut / 2; 
    }
    
    std::cout << "init: " << best_cut << "\n";

    return best_cut;
}


} // namespace pmondriaan
