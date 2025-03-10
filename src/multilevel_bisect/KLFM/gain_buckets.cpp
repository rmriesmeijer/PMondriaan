#include <iostream>
#include <limits>
#include <random>
#include <unordered_set>
#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/KLFM/gain_buckets.hpp"
#include "util/interval.hpp"

namespace pmondriaan {

void gain_buckets::insert(long id, long gain) {
    auto index = gain_to_index(gain);
    if (buckets[index].empty()) {
        if (gain < min_value_present) {
            min_value_present = gain;
        }
        if (gain_to_index(gain) > max_index_present) {
            max_index_present = gain_to_index(gain);
        }
    }
    buckets[index].insert(id);
}

bool gain_buckets::remove(long id, long gain) {
    auto index = gain_to_index(gain);
    if (buckets[index].erase(id) > 0) {
        return true;
    } else {
        return false;
    }
}

long gain_buckets::next() {
    auto index = find_next_index();
    if (index != -1) {
        return *buckets[index].begin();
    } else {
        return -1;
    }
}

long gain_buckets::gain_next() {
    auto index = find_next_index();
    if (index != -1) {
        return index_to_gain(index);
    } else {
        return std::numeric_limits<long>::min();
    }
}

void gain_buckets::print() {
    std::cout << "Buckets: \n";
    for (auto i = 0u; i < buckets.size(); i++) {
        std::cout << index_to_gain(i) << ": ";
        for (auto j : buckets[i]) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }
}

long gain_buckets::find_next_index() {
    auto index = max_index_present;
    if (index == -1) {
        return -1;
    }
    long min_index = gain_to_index(min_value_present);
    while ((index >= min_index) && buckets[index].empty()) {
        index--;
    }
    max_index_present = index;
    if (index >= min_index) {
        return index;
    } else {
        return -1;
    }
}

void gain_structure::init_() {
    auto max_gain = std::vector<long>(2, buckets[0].index_to_gain(0));
    for (auto i = 0u; i < H_.size(); i++) {
        auto& v = H_(i);
        long gain = 0;
        long from = v.part();
        long to = (v.part() + 1) % 2;
        for (auto n : v.nets()) {
            if (C_[H_.local_id_net(n)][from] == 1) {
                gain += H_.net(n).cost();
            }
            if (C_[H_.local_id_net(n)][to] == 0) {
                gain -= H_.net(n).cost();
            }
        }
        buckets[from].insert(v.id(), gain);
        gains[i] = gain;
        if (gain > max_gain[from]) {
            max_gain[from] = gain;
        }
    }

    buckets[0].set_max_present(buckets[0].gain_to_index(max_gain[0]));
    buckets[1].set_max_present(buckets[1].gain_to_index(max_gain[1]));
}

long gain_structure::part_next(long max_extra_weight_0, long max_extra_weight_1, std::mt19937& rng) {
    auto v0 = buckets[0].next();
    auto v1 = buckets[1].next();
    if (v0 == -1) {
        return 1;
    }
    if (v1 == -1) {
        return 0;
    }
    auto gain_v0 = std::numeric_limits<long>::min();
    auto gain_v1 = std::numeric_limits<long>::min();
    if (max_extra_weight_0 - H_(H_.local_id(v0)).weight() > 0) {
        gain_v0 = buckets[0].gain_next();
    }
    if (max_extra_weight_1 - H_(H_.local_id(v1)).weight() > 0) {
        gain_v1 = buckets[1].gain_next();
    }
    if (gain_v0 > gain_v1) {
        return 0;
    }
    if (gain_v1 > gain_v0) {
        return 1;
    }
    return rng() % 2;
}

void gain_structure::move(long v) {
    auto& vertex = H_(H_.local_id(v));
    long from = vertex.part();
    long to = (vertex.part() + 1) % 2;

    // We move v and remove v from its bucket
    H_.move_sorted(v, C_);
    buckets[from].remove(v, gains[H_.local_id(v)]);
    gains[H_.local_id(v)] = std::numeric_limits<long>::min();

    for (auto n : vertex.nets()) {
        if (C_[H_.local_id_net(n)][to] == 0) {
            for (auto u : H_.net(n).vertices()) {
                add_gain(u, H_.net(n).cost());
            }
        }
        if (C_[H_.local_id_net(n)][to] == 1) {
            long u = -1;
            if (to == 1) {
                u = H_.net(n).vertices()[H_.net(n).size() - 1];
            } else {
                u = H_.net(n).vertices()[0];
            }

            if (H_(H_.local_id(u)).part() != to || u == v) {
                std::cout << "Adding gain to wrong vertex!!\n";
                for (auto t : H_.net(n).vertices()) {
                    std::cout << H_(H_.local_id(t)).part() << " ";
                }
            }
            add_gain(u, -1 * H_.net(n).cost());
        }

        C_[H_.local_id_net(n)][to]++;
        C_[H_.local_id_net(n)][from]--;

        if (C_[H_.local_id_net(n)][from] == 0) {
            for (auto u : H_.net(n).vertices()) {
                add_gain(u, -1 * H_.net(n).cost());
            }
        }
        if (C_[H_.local_id_net(n)][from] == 1) {
            long u = -1;
            if (from == 1) {
                u = H_.net(n).vertices()[H_.net(n).size() - 1];
            } else {
                u = H_.net(n).vertices()[0];
            }

            if (H_(H_.local_id(u)).part() != from || u == v) {
                std::cout << "Adding gain to wrong vertex from " << from
                          << "C: " << C_[H_.local_id_net(n)][0] << " "
                          << C_[H_.local_id_net(n)][1] << "!!\n";
                for (auto t : H_.net(n).vertices()) {
                    std::cout << t << " " << H_(H_.local_id(t)).part() << " ";
                }
            }
            add_gain(u, H_.net(n).cost());
        }
    }
}

void gain_structure::move(long v, std::vector<std::vector<long>>& C_loc) {
    auto& vertex = H_(H_.local_id(v));
    long from = vertex.part();
    long to = (vertex.part() + 1) % 2;

    // We move v and remove v from its bucket
    H_.move_sorted(v, C_loc);
    buckets[from].remove(v, gains[H_.local_id(v)]);
    gains[H_.local_id(v)] = std::numeric_limits<long>::min();

    for (auto n : vertex.nets()) {
        if (C_[H_.local_id_net(n)][to] == 0) {
            for (auto u : H_.net(n).vertices()) {
                add_gain(u, H_.net(n).cost());
            }
        }
        if ((C_[H_.local_id_net(n)][to] == 1) && (C_loc[H_.local_id_net(n)][to] == 1)) {
            long u = -1;
            if (to == 1) {
                u = H_.net(n).vertices()[H_.net(n).size() - 1];
            } else {
                u = H_.net(n).vertices()[0];
            }

            if (H_(H_.local_id(u)).part() != to || u == v) {
                std::cout << "Adding gain to wrong vertex!!";
            }
            add_gain(u, -1 * H_.net(n).cost());
        }

        C_[H_.local_id_net(n)][to]++;
        C_[H_.local_id_net(n)][from]--;
        C_loc[H_.local_id_net(n)][to]++;
        C_loc[H_.local_id_net(n)][from]--;

        if (C_[H_.local_id_net(n)][from] == 0) {
            for (auto u : H_.net(n).vertices()) {
                add_gain(u, -1 * H_.net(n).cost());
            }
        }
        if ((C_[H_.local_id_net(n)][from] == 1) && (C_loc[H_.local_id_net(n)][from] == 1)) {
            long u = -1;
            if (from == 1) {
                u = H_.net(n).vertices()[H_.net(n).size() - 1];
            } else {
                u = H_.net(n).vertices()[0];
            }

            if (H_(H_.local_id(u)).part() != from || u == v) {
                std::cout << "Adding gain to wrong vertex!!";
            }
            add_gain(u, H_.net(n).cost());
        }
    }
}

void gain_structure::remove(long v) {
    auto& vertex = H_(H_.local_id(v));
    long from = vertex.part();
    if (!buckets[from].remove(v, gains[H_.local_id(v)])) {
        std::cerr << "Error: Could not remove v from buckets";
    }
    gains[H_.local_id(v)] = std::numeric_limits<long>::min();
}

bool gain_structure::done() {
    if ((buckets[0].next() == -1) && (buckets[1].next() == -1)) {
        return true;
    } else {
        return false;
    }
}


long gain_structure::compute_size_buckets() {
    long max_value = 0;
    for (auto& v : H_.vertices()) {
        long sum = 0;
        for (auto n : v.nets()) {
            sum += H_.net(n).cost();
        }
        if (sum > max_value) {
            max_value = sum;
        }
    }
    return 2 * max_value + 1;
}

void gain_structure::add_gain(long v, long value) {
    long local_id = H_.local_id(v);
    long old_gain = gains[local_id];
    if (old_gain == std::numeric_limits<long>::min()) {
        return;
    }
    buckets[H_(local_id).part()].remove(v, old_gain);
    long new_gain = old_gain + value;
    buckets[H_(local_id).part()].insert(v, new_gain);
    gains[local_id] = new_gain;
}

void gain_structure::check_gains() {
    for (auto i = 0u; i < H_.size(); i++) {
        if (gains[i] != std::numeric_limits<long>::min()) {
            auto& v = H_(i);
            long gain = 0;
            long from = v.part();
            long to = (v.part() + 1) % 2;
            for (auto n : v.nets()) {
                if (C_[H_.local_id_net(n)][from] == 1) {
                    gain += H_.net(n).cost();
                }
                if (C_[H_.local_id_net(n)][to] == 0) {
                    gain -= H_.net(n).cost();
                }
            }
            if (gain != gains[i]) {
                std::cout << "Gain of vertex " << v.id() << " incorrect! Is "
                          << gains[i] << " should be " << gain << "\n";
            }
        }
    }
}

} // namespace pmondriaan
