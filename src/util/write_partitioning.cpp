#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include "util/write_partitioning.hpp"

namespace pmondriaan {

bool partitioning_to_file(bulk::world& world, pmondriaan::hypergraph& H, std::string file, int k) {
    auto starts = start_parts(world, H, k);
    if (world.rank() == 0) {
        std::ofstream out2("./fixfile_x.txt");
        if (out2.fail()) {
            std::cerr << "Error: " << std::strerror(errno);
            return false;
        }
        std::ofstream out4("./fixfile_o.txt");
        if (out4.fail()) {
            std::cerr << "Error: " << std::strerror(errno);
            return false;
        }
        for(auto v : H.vertices()) {
            out2 << v.id() << " " << v.part() << "\n";
            if(rand() % 1000 < 3) {
                out4 << v.id() << " " << v.part() << "\n";
            }
            else {
                out4 << v.id() << " -1\n";
            }
        }
        std::ofstream out3("./fixfile_x2.txt");
        if (out3.fail()) {
            std::cerr << "Error: " << std::strerror(errno);
            return false;
        }
        std::ofstream out5("./fixfile_o2.txt");
        if (out5.fail()) {
            std::cerr << "Error: " << std::strerror(errno);
            return false;
        }
        for(auto v : H.vertices()) {
            out3 << v.part() << "\n";
            if(rand() % 1000 < 3) {
                out5 << v.part() << "\n";
            }
            else {
                out5 << "-1\n";
            }
        }
        out2.close();
        out3.close();
        std::ofstream out(file);
        if (out.fail()) {
            std::cerr << "Error: " << std::strerror(errno);
            return false;
        }
        out << "%%MatrixMarket distributed-matrix coordinate pattern general\n";
        out << H.global_number_nets() << " " << H.global_size() << " "
            << H.nr_nz() << " " << k << "\n";
        for (auto i : starts) {
            out << i << "\n";
        }
        out.close();
    }
    for (int i = 0; i < k; i++) {
        for (long turn = 0; turn < world.active_processors(); turn++) {
            if (turn == world.rank()) {
                std::ofstream out(file, std::ios::app);
                if (out.fail()) {
                    std::cerr << "Error: " << std::strerror(errno);
                    return false;
                }
                for (auto& v : H.vertices()) {
                    if (v.part() == i) {
                        for (auto n : v.nets()) {
                            out << n + 1 << " " << v.id() + 1 << "\n";
                        }
                    }
                }
                out.close();
            }
            world.sync();
        }
    }
    return true;
}

std::vector<long> start_parts(bulk::world& world, pmondriaan::hypergraph& H, int k) {
    auto counts = bulk::coarray<long>(world, k);
    for (int i = 0; i < k; i++) {
        counts[i] = 0;
    }
    for (auto& v : H.vertices()) {
        counts[v.part()] += v.degree();
    }
    auto result =
    bulk::foldl_each(counts, [](auto& lhs, auto rhs) { lhs += rhs; });
    result.insert(result.begin(), 0);
    for (int i = 1; i < k + 1; i++) {
        result[i + 1] += result[i];
    }
    return result;
}

} // namespace pmondriaan