#include <array>
#include <limits>
#include <queue>
#include <random>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/KLFM/KLFM_parallel.hpp"
#include "multilevel_bisect/KLFM/gain_buckets.hpp"

namespace pmondriaan {

/**
 * Runs the KLFM algorithm in parallel to improve a given partitioning. Return the cutsize of the best solution.
 */
long KLFM_par(bulk::world& world,
              pmondriaan::hypergraph& H,
              std::vector<std::vector<long>>& C,
              long weight_0,
              long weight_1,
              long max_weight_0,
              long max_weight_1,
              pmondriaan::options& opts,
              std::mt19937& rng,
              long cut_size) {

    size_t pass = 0;
    long prev_cut_size;

    // TODO: use C to compute cutsize without communication
    if (cut_size == std::numeric_limits<decltype(cut_size)>::max()) {
        prev_cut_size = pmondriaan::cutsize(world, H, opts.metric);
    } else {
        prev_cut_size = cut_size;
    }
    auto total_weights = std::array<long, 2>();
    total_weights[0] = weight_0;
    total_weights[1] = weight_1;

    while (pass < opts.KLFM_max_passes) {
        auto result = KLFM_pass_par(world, H, C, prev_cut_size, total_weights,
                                    max_weight_0, max_weight_1, opts, rng);
        if (result < prev_cut_size) {
            prev_cut_size = result;
        } else {
            break;
        }
        pass++;
    }
    return prev_cut_size;
}

/**
 * Runs a single pass of the KLFM algorithm in parallel to improve a given partitioning.
 */
long KLFM_pass_par(bulk::world& world,
                   pmondriaan::hypergraph& H,
                   std::vector<std::vector<long>>& C,
                   long cut_size,
                   std::array<long, 2>& total_weights,
                   long max_weight_0,
                   long max_weight_1,
                   pmondriaan::options& opts,
                   std::mt19937& rng) {
    auto s = world.rank();
    auto p = world.active_processors();

    auto gain_structure = pmondriaan::gain_structure(H, C);

    long best_cut_size = cut_size;
    auto no_improvement_moves = std::vector<long>();

    // Each processor is responsible for keeping track of some nets
    auto net_partition =
    bulk::block_partitioning<1>({H.global_number_nets()}, {(size_t)p});
    /*We keep track of the previous counts we are responsible for,
    we have the count of part 0 and the count of part 1 for each net */
    auto previous_C = bulk::coarray<long>(world, net_partition.local_count(s) * 2);
    auto cost_my_nets = bulk::coarray<long>(world, net_partition.local_count(s));
    auto cut_size_my_nets =
    init_previous_C(world, H, C, previous_C, cost_my_nets, net_partition);

    // Stores the proposed moves as: gain, weight change, processor id
    auto moves_queue = bulk::queue<long, long, int, long>(world);
    auto update_nets = bulk::queue<long, long>(world);
    bulk::var<long> rejected(world);
    bulk::var<long> new_weight_0(world);
    bulk::var<long> new_weight_1(world);
    bulk::var<int> done(world);
    done = 0;

    bool all_done = false;
    while (!all_done) {
        auto prev_total_weights = total_weights;

        // Find best KLFM_par_number_send_moves moves
        // A move is stored as v_to_move, gain_v, weight change
        auto moves =
        std::vector<std::tuple<long, long, long>>(opts.KLFM_par_number_send_moves,
                                                  std::make_tuple(-1, 0, 0));
        find_top_moves(H, gain_structure, moves, total_weights, max_weight_0,
                       max_weight_1, rng);

        // Send moves to processor 0
        long tag = 0;
        for (auto& move : moves) {
            if (std::get<0>(move) == -1) {
                break;
            }
            moves_queue(0).send(std::get<1>(move), std::get<2>(move), s, tag);
            tag++;
        }
        world.sync();

        if (s == 0) {
            auto rejected_all = reject_unbalanced_moves(p, moves_queue, prev_total_weights,
                                                        max_weight_0, max_weight_1);
            for (auto t = 0; t < p; t++) {
                rejected(t) = rejected_all[t];
            }
            new_weight_0.broadcast(prev_total_weights[0]);
            new_weight_1.broadcast(prev_total_weights[1]);
        }
        world.sync();
        total_weights[0] = new_weight_0.value();
        total_weights[1] = new_weight_1.value();

        // We reverse moves that have not been selected
        auto index = moves.size() - 1;
        if (rejected.value() > 0) {
            while (rejected.value() != 0) {
                if (std::get<2>(moves[index]) > 0) {
                    auto& vertex = H(H.local_id(std::get<0>(moves[index])));
                    H.move(vertex.id(), C);
                    rejected--;
                }
                index--;
            }
        } else {
            while (rejected.value() != 0) {
                if (std::get<2>(moves[index]) < 0) {
                    auto& vertex = H(H.local_id(std::get<0>(moves[index])));
                    H.move(vertex.id(), C);
                    rejected++;
                }
                index--;
            }
        }

        // We send updates about counts in nets to responsible processors
        // TODO: Make this more efficient by not sending everything
        for (auto i = 0u; i < H.nets().size(); i++) {
            update_nets(net_partition.owner(H.nets()[i].id()))
            .send(H.nets()[i].id(), C[H.local_id_net(H.nets()[i].id())][0]);
        }
        world.sync();

        cut_size_my_nets =
        update_C(world, H, C, previous_C, update_nets, net_partition,
                 cost_my_nets, cut_size_my_nets, true, gain_structure);

        // We also send all processors the total cutsize of the nets this p is responsible for
        cut_size = bulk::sum(world, cut_size_my_nets);

        if (cut_size > best_cut_size) {
            for (auto& move : moves) {
                no_improvement_moves.push_back(std::get<0>(move));
            }
        } else {
            best_cut_size = cut_size;
            no_improvement_moves.clear();
        }

        if (gain_structure.done()) {
            done = 1;
        }
        if (bulk::sum(done) == world.active_processors()) {
            all_done = true;
        }
    }

    // Reverse moves up to best point
    long weight_change = 0;
    auto changed_nets = std::unordered_set<long>();
    for (auto v : no_improvement_moves) {
        auto& vertex = H(H.local_id(v));
        if (vertex.part() == 0) {
            weight_change -= vertex.weight();
        } else {
            weight_change += vertex.weight();
        }
        for (auto n : vertex.nets()) {
            changed_nets.insert(n);
        }
        H.move(v, C);
    }

    // We send updates about counts in nets to responsible processors
    for (auto n : changed_nets) {
        update_nets(net_partition.owner(n)).send(n, C[H.local_id_net(n)][0]);
    }
    world.sync();
    update_C(world, H, C, previous_C, update_nets, net_partition, cost_my_nets,
             cut_size_my_nets, false, gain_structure);

    auto total_change = bulk::sum(world, weight_change);
    total_weights[0] += total_change;
    total_weights[1] -= total_change;

    return best_cut_size;
}

/**
 * Initializes the previous_C counts using communication and return the cutsize
 * of the nets the processor is responsible for.
 */
long init_previous_C(bulk::world& world,
                     pmondriaan::hypergraph& H,
                     std::vector<std::vector<long>>& C,
                     bulk::coarray<long>& previous_C,
                     bulk::coarray<long>& cost_my_nets,
                     bulk::block_partitioning<1>& net_partition) {
    auto s = world.rank();
    for (auto i = 0u; i < net_partition.local_count(s); i++) {
        cost_my_nets[i] = 0;
    }
    for (auto i = 0u; i < net_partition.local_count(s) * 2; i++) {
        previous_C[i] = 0;
    }

    for (auto i = 0u; i < C.size(); i++) {
        auto net = H.global_id_net(i);
        previous_C(net_partition.owner(net))[2 * net_partition.local(net)[0]] = C[i][0];
        previous_C(net_partition.owner(net))[(2 * net_partition.local(net)[0]) + 1] =
        C[i][1];
        cost_my_nets(net_partition.owner(net))[net_partition.local(net)[0]] =
        H.net(net).cost();
    }
    world.sync();

    long cut_size_my_nets = 0;
    for (auto i = 0u; i < net_partition.local_count(s); i++) {
        if (!((previous_C[i * 2] == 0) || (previous_C[(i * 2) + 1] == 0))) {
            cut_size_my_nets += cost_my_nets[i];
        }
    }
    return cut_size_my_nets;
}

/**
 * Updates the gain values that were outdated.
 */
void update_gains(pmondriaan::hypergraph& H,
                  pmondriaan::net& net,
                  std::vector<long> C_loc,
                  std::vector<long> C_new,
                  pmondriaan::gain_structure& gain_structure) {
    if (((C_new[0] == 0) && (C_loc[0] > 0)) || ((C_new[1] == 0) && (C_loc[1] > 0))) {
        for (auto v : net.vertices()) {
            gain_structure.add_gain(v, -1 * net.cost());
        }
    }
    if ((C_new[0] == 1) && (C_loc[0] > 1)) {
        for (auto v : net.vertices()) {
            if (H(H.local_id(v)).part() == 0) {
                gain_structure.add_gain(v, net.cost());
            }
        }
    }
    if ((C_new[1] == 1) && (C_loc[1] > 1)) {
        for (auto v : net.vertices()) {
            if (H(H.local_id(v)).part() == 1) {
                gain_structure.add_gain(v, net.cost());
            }
        }
    }
    if (((C_new[0] > 0) && (C_loc[0] == 0)) || ((C_new[1] > 0) && (C_loc[1] == 0))) {
        for (auto v : net.vertices()) {
            gain_structure.add_gain(v, net.cost());
        }
    }
    if ((C_new[0] > 1) && (C_loc[0] == 1)) {
        for (auto v : net.vertices()) {
            if (H(H.local_id(v)).part() == 0) {
                gain_structure.add_gain(v, -1 * net.cost());
            }
        }
    }
    if ((C_new[1] > 1) && (C_loc[1] == 1)) {
        for (auto v : net.vertices()) {
            if (H(H.local_id(v)).part() == 1) {
                gain_structure.add_gain(v, -1 * net.cost());
            }
        }
    }
}

/**
 * Finds the best moves for a processor sequentially, by only updating local data.
 */
void find_top_moves(pmondriaan::hypergraph& H,
                    pmondriaan::gain_structure& gain_structure,
                    std::vector<std::tuple<long, long, long>>& moves,
                    std::array<long, 2>& weights,
                    long max_weight_0,
                    long max_weight_1,
                    std::mt19937& rng) {
    auto max_extra_weight = std::array<long, 2>();
    auto moves_found = 0u;

    while (!gain_structure.done() && (moves_found < moves.size())) {
        max_extra_weight[0] = max_weight_0 - weights[0];
        max_extra_weight[1] = max_weight_1 - weights[1];

        long part_to_move =
        gain_structure.part_next(max_extra_weight[0], max_extra_weight[1], rng);
        auto v_to_move = gain_structure.next(part_to_move);

        if (max_extra_weight[(part_to_move + 1) % 2] - H(H.local_id(v_to_move)).weight() >= 0) {
            auto weight_v = H(H.local_id(v_to_move)).weight();
            if (part_to_move == 0) {
                moves[moves_found] =
                std::make_tuple(v_to_move, gain_structure.gain_next(part_to_move),
                                -1 * weight_v);
                weights[0] -= weight_v;
                weights[1] += weight_v;
            } else {
                moves[moves_found] =
                std::make_tuple(v_to_move, gain_structure.gain_next(part_to_move), weight_v);
                weights[1] -= weight_v;
                weights[0] += weight_v;
            }
            gain_structure.move(v_to_move);
            moves_found++;
        } else {
            gain_structure.remove(v_to_move);
        }
    }
}

/**
 * Determines how many moves from part 0 or 1 should be rejected on each processor.
 * Return the number of rejected moves for each processor. This is positive when
 * we have to move back vertices from part 0 to part 1 and negative otherwise.
 */
std::vector<long> reject_unbalanced_moves(long p,
                                          bulk::queue<long, long, int, long>& moves_queue,
                                          std::array<long, 2>& total_weights,
                                          long max_weight_0,
                                          long max_weight_1) {
    // The number of moves to be rejected for each processor
    auto rejected = std::vector<long>(p, 0);
    auto done = std::vector<bool>(p, false);

    // Sort queue on gain
    std::sort(moves_queue.begin(), moves_queue.end(), [](auto lhs, auto rhs) {
        if (std::get<0>(lhs) == std::get<0>(rhs)) {
            return std::get<3>(lhs) > std::get<3>(rhs);
        } else {
            return std::get<0>(lhs) < std::get<0>(rhs);
        }
    });
    std::queue<std::pair<long, long>> received_moves;
    long total_balance = 0;
    for (const auto& [gain, weight_change, t, tag] : moves_queue) {
        total_balance += weight_change;
        received_moves.push(std::make_pair(weight_change, t));
    }
    total_weights[0] += total_balance;
    total_weights[1] -= total_balance;
    while (total_weights[0] > max_weight_0) {
        if (received_moves.empty()) {
            break;
        }
        auto weight_change = received_moves.front().first;
        auto t = received_moves.front().second;
        if ((weight_change > 0) && (!done[t])) {
            if (total_weights[1] + weight_change <= max_weight_1) {
                rejected[t]++;
                total_weights[0] -= weight_change;
                total_weights[1] += weight_change;
            } else {
                done[t] = true;
            }
        }
        received_moves.pop();
    }

    while (total_weights[1] > max_weight_1) {
        if (received_moves.empty()) {
            break;
        }
        auto weight_change = received_moves.front().first;
        auto t = received_moves.front().second;
        if ((weight_change < 0) && (!done[t])) {
            if (total_weights[0] - weight_change <= max_weight_0) {
                rejected[t]--;
                total_weights[0] -= weight_change;
                total_weights[1] += weight_change;
            } else {
                done[t] = true;
            }
        }
        received_moves.pop();
    }
    return rejected;
}

/**
 * Updates the counts by communicating with the responsible processor, returns the new cutsize_my_nets.
 */
long update_C(bulk::world& world,
              pmondriaan::hypergraph& H,
              std::vector<std::vector<long>>& C,
              bulk::coarray<long>& previous_C,
              bulk::queue<long, long>& update_nets,
              bulk::partitioning<1>& net_partition,
              bulk::coarray<long>& cost_my_nets,
              long cut_size_my_nets,
              bool update_g,
              pmondriaan::gain_structure& gain_structure) {
    auto s = world.rank();
    auto C_new = std::vector<std::vector<long>>(net_partition.local_count(s),
                                                std::vector<long>(2, 0));
    for (auto i = 0u; i < net_partition.local_count(s); i++) {
        C_new[i][0] = previous_C[2 * i];
        C_new[i][1] = previous_C[(2 * i) + 1];
    }
    // Update the counts of my nets
    for (const auto& [net, C_0] : update_nets) {
        auto local = net_partition.local(net)[0];
        auto change = C_0 - previous_C[2 * local];
        if (change != 0) {
            if ((C_new[local][0] == 0) || (C_new[local][1] == 0)) {
                cut_size_my_nets += cost_my_nets[net_partition.local(net)[0]];
            }
            C_new[local][0] += change;
            C_new[local][1] -= change;
            // TODO make sure we also adjust the gains here
            if ((C_new[local][0] == 0) || (C_new[local][1] == 0)) {
                cut_size_my_nets -= cost_my_nets[net_partition.local(net)[0]];
            }
        }
    }
    for (auto i = 0u; i < net_partition.local_count(s); i++) {
        previous_C[2 * i] = C_new[i][0];
        previous_C[(2 * i) + 1] = C_new[i][1];
    }

    // Each processor uses a get to find the updated counts it needs.
    auto future_counts = std::vector<bulk::future<long>>();
    for (auto i = 0u; i < H.nets().size(); i++) {
        auto net = H.nets()[i].id();
        future_counts.push_back(
        previous_C(net_partition.owner(net))[2 * net_partition.local(net)[0]].get());
    }
    world.sync();

    if (update_g) {
        for (auto i = 0u; i < H.nets().size(); i++) {
            if (future_counts[i].value() != C[i][0]) {
                update_gains(H, H.nets()[i], C[i],
                             std::vector<long>({future_counts[i].value(),
                                                (long)H.nets()[i].global_size() -
                                                future_counts[i].value()}),
                             gain_structure);
            }
        }
    }
    for (auto i = 0u; i < H.nets().size(); i++) {
        C[i][0] = future_counts[i].value();
        C[i][1] = H.nets()[i].global_size() - future_counts[i].value();
    }
    return cut_size_my_nets;
}

// For testing purposes
void check_C(bulk::world& world, pmondriaan::hypergraph& H, std::vector<std::vector<long>>& C) {
    auto correct_C = init_counts(world, H);
    for (auto i = 0u; i < C.size(); i++) {
        if (C[i][0] != correct_C[i][0]) {
            world.log("s %d: C[%d][0] incorrect", world.rank(), i);
        }
        if (C[i][1] != correct_C[i][1]) {
            world.log("s %d: C[%d][1] incorrect", world.rank(), i);
        }
    }
}

} // namespace pmondriaan
