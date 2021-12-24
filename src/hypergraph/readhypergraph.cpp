#include "hypergraph/readhypergraph.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include <hypergraph/hypergraph.hpp>

namespace pmondriaan {

/**
 * Creates a hypergraph from a graph in mtx format.
 */
std::optional<pmondriaan::hypergraph>
read_hypergraph_istream(std::istream& fin, std::istream& ffin, std::string mode_weight) {

    size_t E, V;
    uint64_t L;
    uint64_t nz = 0;
    std::string line;
    std::getline(fin, line);
    std::string object, format, field, symmetry;
    std::istringstream iss(line);
    if (!(iss >> object >> object >> format >> field >> symmetry)) {
        return std::nullopt;
    }

    // Ignore headers and comments:
    while (fin.peek() == '%') {
        fin.ignore(2048, '\n');
    }

    // Read defining parameters:
    fin >> E >> V >> L;

    auto nets_list = std::vector<std::vector<long>>(V);
    auto vertex_list = std::vector<std::vector<long>>(E);

    // Read the data
    std::getline(fin, line);

    if (symmetry == "general") {
        nz = L;
        while (std::getline(fin, line)) {
            size_t e, v;
            std::istringstream iss(line);
            if (!(iss >> e >> v)) {
                return std::nullopt;
            }
            nets_list[v - 1].push_back(e - 1);
            vertex_list[e - 1].push_back(v - 1);
        }
    } else {
        if (symmetry == "symmetric") {
            while (std::getline(fin, line)) {
                size_t e, v;
                nz++;
                std::istringstream iss(line);
                if (!(iss >> e >> v)) {
                    return std::nullopt;
                }
                nets_list[v - 1].push_back(e - 1);
                vertex_list[e - 1].push_back(v - 1);
                if (v != e) {
                    nz++;
                    nets_list[e - 1].push_back(v - 1);
                    vertex_list[v - 1].push_back(e - 1);
                }
            }
        } else {
            std::cerr << "Error: unknown symmetry";
            return std::nullopt;
        }
    }

    auto vertices = std::vector<pmondriaan::vertex>();
    auto nets = std::vector<pmondriaan::net>();
    if (mode_weight == "one") {
        for (size_t i = 0; i < V; i++) {
            vertices.push_back(pmondriaan::vertex(i, nets_list[i]));
        }
    } else if (mode_weight == "degree") {
        for (size_t i = 0; i < V; i++) {
            vertices.push_back(pmondriaan::vertex(i, nets_list[i], nets_list[i].size()));
        }
    } else {
        std::cerr << "Error: unknown mode_weight";
        return std::nullopt;
    }

    for (size_t i = 0; i < E; i++) {
        if (!vertex_list[i].size()) {
            nets.push_back(pmondriaan::net(i, vertex_list[i]));
        }
    }


    auto H = pmondriaan::hypergraph(V, E, vertices, nets, nz);

    // List of fixed vertices of each bisection part.
    auto fixed1 = std::vector<long>();
    auto fixed2 = std::vector<long>();
    std::string ffline;

    // Read fixed file
    while (std::getline(ffin, ffline)) {
        try {
            std::vector<long> pair = std::vector<long>();
            char *pch;
            const char s[2] = " ";
            pch = strtok ((char*)ffline.c_str(),s);
            while (pch != NULL) {
                pair.push_back(std::stol(pch));
                pch = strtok (NULL, " ");
            }
            if(pair[1] == 0) {
                fixed1.push_back(pair[0]);
                H(H.local_id(pair[0])).set_fixpart(0);
            }
            if(pair[1] ==  1){
                fixed2.push_back(pair[0]);
                H(H.local_id(pair[0])).set_fixpart(1);
            }
        } catch (...) { std::cerr << "Error: fixfile\n"; }
    }

    // We need to transform due to fixed vertices in both parts.
    if(fixed1.size() != 0 && fixed2.size() != 0) {
        std::cout << "transform fixed\n";
        auto fix_nets1 = std::vector<long>();
        auto fix_nets2 = std::vector<long>();
        auto fix_nets_overlap = std::vector<long>();
        auto netsf = std::vector<long>();
        long w1 = 0;
        long idf = fixed1[0];
        long g = 0;
        long maxid = 0;

        long fc = 0;
        for(auto& n : H.nets()) {
            auto vec = std::vector<long>();
            bool net_has_fix1 = false;
            bool net_has_fix2 = false;
            if(n.id() > maxid) {
                maxid = n.id();
            }
            for(auto v : n.vertices()) {
                if(H(H.local_id(v)).fixpart() == 0) {
                    net_has_fix1 = true;
                } else if(H(H.local_id(v)).fixpart() == 1) {
                    net_has_fix2 = true;
                } else {
                    vec.push_back(v);
                }
            }
            if(net_has_fix1 == true && net_has_fix2 == true) {
                fc += n.cost();
                fix_nets_overlap.push_back(n.id());
            } else if(net_has_fix2 == true) {
                fix_nets2.push_back(n.id());
                netsf.push_back(n.id());
            } else if(net_has_fix1 == true) {
                fix_nets1.push_back(n.id());
                netsf.push_back(n.id());
            }

            while(n.vertices().size() > 0) {
                n.pop_back();
            }
            for(auto indi : vec) {
                n.add_vertex(indi);
            }
            
            g++;
        }
        std::cout << "transform fixed1\n";

        for(auto n : fix_nets_overlap) {
            H.remove_net_by_index(H.local_id_net(n));
        }

        for(auto i : fixed1) {
            w1 += H(H.local_id(i)).weight();
            H.remove_free_vertex(i);
            g++;
        }
        std::cout << "transform fixed2\n";

        long w2 = 0;
        for(auto i : fixed2) {
            w2 += H(H.local_id(i)).weight();
            H.remove_free_vertex(i);
            g++;
        }

        std::cout << "transform fixed3 " << H.nets().size() << "\n";

        long fix_weight_difference = 0;
        long fixed_vertex_part = -1;
        long wf = -1;
        long cid = maxid;
        if(w1 > w2) {
            wf = w1 - w2;
            fix_weight_difference = w2;
            fixed_vertex_part = 0;
            for(auto n : fix_nets2) {
                if(H.nets()[H.local_id_net(n)].size() > 0) {
                    long costc = H.nets()[H.local_id_net(n)].cost();
                    auto verts = std::vector<long>();
                    for(auto i : H.nets()[H.local_id_net(n)].vertices()) {
                        verts.push_back(i);
                    }
                    if(verts.size() == 0) {
                        std::cerr << "bad1! \n"; 
                    }
                    H.nets()[H.local_id_net(n)].set_cost( costc);
                    //H.add_net(cid, verts, costc);
                    cid++;
                }
            }
        }
        else {
            wf = w2 - w1;
            fix_weight_difference = w1;
            fixed_vertex_part = 1;
            for(auto n : fix_nets1) {
                if(H.nets()[H.local_id_net(n)].size() > 0) {
                    long costc = H.nets()[H.local_id_net(n)].cost();
                    auto verts = std::vector<long>();
                    for(auto i : H.nets()[H.local_id_net(n)].vertices()) {
                        verts.push_back(i);
                    }
                    if(verts.size() == 0) {
                        std::cerr << "bad! \n"; 
                    }
                    H.nets()[H.local_id_net(n)].set_cost( costc);
                    //H.add_net(cid, verts, costc);
                    cid++;
                }
            }
        }

        H.add_vertex(idf, netsf, wf);
        H.add_to_nets(H(H.local_id(idf)));

        auto val = 0;
        for(auto n : H.nets()) {
            val += n.size();
        }
        size_t nzx = val;
        size_t VX = V - fixed1.size() - fixed2.size() + 1;
        size_t EX = H.nets().size();
        auto verticesx = std::vector<pmondriaan::vertex>();
        auto netsx = std::vector<pmondriaan::net>();



        for (size_t i = 0; i < H.vertices().size(); i++) {
            for(size_t j = 0; j < H.vertices()[i].nets().size(); j++) {
                H.vertices()[i].nets()[j] = H.local_id_net(H.vertices()[i].nets()[j]);
            }
            verticesx.push_back(pmondriaan::vertex(i, H.vertices()[i].nets(), H.vertices()[i].weight()));
        }

        for(auto& n : H.nets()) {
            for(size_t i = 0; i < n.vertices().size(); i++) {
                n.vertices()[i] = H.local_id(n.vertices()[i]);
            }
        }

        for (size_t i = 0; i < EX; i++) {
            if (H.nets()[i].size() > 0) {
                netsx.push_back(pmondriaan::net(i, H.nets()[i].vertices(), H.nets()[i].cost()));
            }
        }

        auto HC = pmondriaan::hypergraph(VX, EX, verticesx, netsx, nzx);

        HC(H.local_id(idf)).set_fixpart(fixed_vertex_part);
        HC(H.local_id(idf)).set_fixwdiff(fix_weight_difference);
        HC(H.local_id(idf)).set_fixcdiff(fc);

        std::cout << "lol: " << fc << "\n";


        pmondriaan::remove_free_nets(HC, 0);


        
        std::cout << "transform over\n";
        return std::move(HC);
    }


    pmondriaan::remove_free_nets(H, 0);

    std::cout << "transform over 2\n";

    return std::move(H);
}

/**
 * Creates a distributed hypergraph from a graph in mtx format.
 */
std::optional<pmondriaan::hypergraph>
read_hypergraph_istream(std::istream& fin, bulk::world& world, std::istream& ffin, std::string mode_weight) {

    auto s = world.rank();
    auto p = world.active_processors();

    size_t E, V;
    uint64_t L;
    uint64_t nz = 0;

    std::string line;
    std::getline(fin, line);
    std::string object, format, field, symmetry;
    std::istringstream iss(line);
    if (!(iss >> object >> object >> format >> field >> symmetry)) {
        return std::nullopt;
    }

    // Ignore headers and comments:
    while (fin.peek() == '%') {
        fin.ignore(2048, '\n');
    }

    // Read defining parameters:
    fin >> E >> V >> L;

    auto partitioning = bulk::block_partitioning<1>({V}, {(size_t)p});

    // List of nets for each vertex
    auto nets_list = std::vector<std::vector<long>>(partitioning.local_count(s));
    auto vertex_list = std::vector<std::vector<long>>(E);

    // Read the data
    std::getline(fin, line);
    if (symmetry == "general") {
        for (auto i = 0u; i < L; i++) {
            nz = L;
            std::getline(fin, line);
            size_t e, v;
            std::istringstream iss(line);
            if (!(iss >> e >> v)) {
                return std::nullopt;
            } // error
            if (partitioning.owner({v - 1}) == s) {
                long v_loc = partitioning.local({v - 1})[0];
                nets_list[v_loc].push_back(e - 1);
                vertex_list[e - 1].push_back(v - 1);
            }
        }
    } else {
        if (symmetry == "symmetric") {
            for (auto i = 0u; i < L; i++) {
                nz++;
                std::getline(fin, line);
                size_t e, v;
                std::istringstream iss(line);
                if (!(iss >> e >> v)) {
                    return std::nullopt;
                } // error
                if (partitioning.owner({v - 1}) == s) {
                    long v_loc = partitioning.local({v - 1})[0];
                    nets_list[v_loc].push_back(e - 1);
                    vertex_list[e - 1].push_back(v - 1);
                }
                if ((partitioning.owner({e - 1}) == s) && (v != e)) {
                    long v_loc = partitioning.local({e - 1})[0];
                    nets_list[v_loc].push_back(v - 1);
                    vertex_list[v - 1].push_back(e - 1);
                }
                if (v != e) {
                    nz++;
                }
            }
        } else {
            std::cerr << "Error: unknown symmetry";
            return std::nullopt;
        }
    }
    world.sync();

    auto vertices = std::vector<pmondriaan::vertex>();
    auto nets = std::vector<pmondriaan::net>();
    if (mode_weight == "one") {
        for (size_t i = 0; i < partitioning.local_count(s); i++) {
            vertices.push_back(
            pmondriaan::vertex(partitioning.global({i}, s)[0], nets_list[i]));
        }
    } else if (mode_weight == "degree") {
        for (size_t i = 0; i < partitioning.local_count(s); i++) {
            vertices.push_back(pmondriaan::vertex(partitioning.global({i}, s)[0],
                                                  nets_list[i]));
        }
    } else {
        std::cerr << "Error: unknown mode_weight";
        return std::nullopt;
    }
    for (size_t i = 0; i < E; i++) {
        if (!vertex_list[i].empty()) {
            nets.push_back(pmondriaan::net(i, vertex_list[i]));
        }
    }

    auto H = pmondriaan::hypergraph(V, E, vertices, nets, nz);

    if (mode_weight == "degree") {
        // List of fixed vertices of each bisection part.
        auto fixed1 = std::vector<long>();
        auto fixed2 = std::vector<long>();
        std::string ffline;

        // Read fixed file
        while (std::getline(ffin, ffline)) {
            try {
                std::vector<long> pair = std::vector<long>();
                char *pch;
                const char s[2] = " ";
                pch = strtok ((char*)ffline.c_str(),s);
                while (pch != NULL) {
                    pair.push_back(std::stol(pch));
                    pch = strtok (NULL, " ");
                }
                if(pair[1] == 0) {
                    fixed1.push_back(pair[0]);
                    H(H.local_id(pair[0])).set_fixpart(0);
                }
                if(pair[1] ==  1){
                    fixed2.push_back(pair[0]);
                    H(H.local_id(pair[0])).set_fixpart(1);
                }
            } catch (...) { std::cerr << "Error: fixfile\n"; }
        }

        // We need to transform due to fixed vertices in both parts.
        if(fixed1.size() != 0 && fixed2.size() != 0) {
            std::cout << "transform fixed\n";
            auto fix_nets1 = std::vector<long>();
            auto fix_nets2 = std::vector<long>();
            auto fix_nets_overlap = std::vector<long>();
            auto netsf = std::vector<long>();
            long w1 = 0;
            long idf = fixed1[0];
            long g = 0;
            long maxid = 0;

            long fc = 0;
            long fc1 = 0;
            long fc2 = 0;
            for(auto& n : H.nets()) {
                auto vec = std::vector<long>();
                bool net_has_fix1 = false;
                bool net_has_fix2 = false;
                if(n.id() > maxid) {
                    maxid = n.id();
                }
                for(auto v : n.vertices()) {
                    if(H(H.local_id(v)).fixpart() == 0) {
                        net_has_fix1 = true;
                    } else if(H(H.local_id(v)).fixpart() == 1) {
                        net_has_fix2 = true;
                    } else {
                        vec.push_back(v);
                    }
                }
                if(net_has_fix1 == true && net_has_fix2 == true) {
                    fc += n.cost();
                    fix_nets_overlap.push_back(n.id());
                } else if(net_has_fix2 == true) {
                    fc2 += n.cost();
                    fix_nets2.push_back(n.id());
                    netsf.push_back(n.id());
                } else if(net_has_fix1 == true) {
                    fc1 += n.cost();
                    fix_nets1.push_back(n.id());
                    netsf.push_back(n.id());
                }

                while(n.vertices().size() > 0) {
                    n.pop_back();
                }
                for(auto indi : vec) {
                    n.add_vertex(indi);
                }
                
                g++;
            }
            std::cout << "transform fixed1\n";

            for(auto n : fix_nets_overlap) {
                H.remove_net_by_index(H.local_id_net(n));
            }

            for(auto i : fixed1) {
                w1 += H(H.local_id(i)).weight();
                H.remove_free_vertex(i);
                g++;
            }
            std::cout << "transform fixed2\n";

            long w2 = 0;
            for(auto i : fixed2) {
                w2 += H(H.local_id(i)).weight();
                H.remove_free_vertex(i);
                g++;
            }

            std::cout << "transform fixed3 " << H.nets().size() << "\n";

            long fix_weight_difference = 0;
            long fixed_vertex_part = -1;
            long wf = -1;
            long cid = maxid;
            if(w1 > w2) {
                wf = w1 - w2;
                fix_weight_difference = w2;
                fixed_vertex_part = 0;
                fc += fc2;
                for(auto n : fix_nets2) {
                    if(H.nets()[H.local_id_net(n)].size() > 0) {
                        long costc = H.nets()[H.local_id_net(n)].cost();
                        auto verts = std::vector<long>();
                        for(auto i : H.nets()[H.local_id_net(n)].vertices()) {
                            verts.push_back(i);
                        }
                        if(verts.size() == 0) {
                            std::cerr << "bad1! \n"; 
                        }
                        H.nets()[H.local_id_net(n)].set_cost(-1 * costc);
                        H.add_net(cid, verts, costc);
                        cid++;
                    }
                    else {
                        fc--;
                    }
                }
            }
            else {
                wf = w2 - w1;
                fix_weight_difference = w1;
                fixed_vertex_part = 1;
                fc += fc1;
                for(auto n : fix_nets1) {
                    if(H.nets()[H.local_id_net(n)].size() > 0) {
                        long costc = H.nets()[H.local_id_net(n)].cost();
                        auto verts = std::vector<long>();
                        for(auto i : H.nets()[H.local_id_net(n)].vertices()) {
                            verts.push_back(i);
                        }
                        if(verts.size() == 0) {
                            std::cerr << "bad! \n"; 
                        }
                        H.nets()[H.local_id_net(n)].set_cost(-1 * costc);
                        H.add_net(cid, verts, costc);
                        cid++;
                    }
                    else {
                        fc--;
                    }
                }
            }

            long mincost = 0;
            long maxcost = 0;
            for(auto& n : H.nets()) {
                if(mincost > n.cost()) {
                    mincost = n.cost();
                }
                if(maxcost < n.cost()) {
                    maxcost = n.cost();
                }
            }
            std::cout << "transform fixed4 minco " << mincost << "\n";
            std::cout << "transform fixed4 maxco " << maxcost << "\n";

            std::cout << "transform fixed4 " << H.nets().size() << "\n";

            std::cout << "Fix nets: " << fix_nets1.size() + fix_nets2.size() << "\n";

            H.add_vertex(idf, netsf, wf);
            H.add_to_nets(H(H.local_id(idf)));

            auto val = 0;
            for(auto n : H.nets()) {
                val += n.size();
            }

            size_t nzx = val;
            size_t VX = V - fixed1.size() - fixed2.size() + 1;
            size_t EX = H.nets().size();
            auto verticesx = std::vector<pmondriaan::vertex>();
            auto netsx = std::vector<pmondriaan::net>();

            for (size_t i = 0; i < H.vertices().size(); i++) {
                for(size_t j = 0; j < H.vertices()[i].nets().size(); j++) {
                    H.vertices()[i].nets()[j] = H.local_id_net(H.vertices()[i].nets()[j]);
                }
                verticesx.push_back(pmondriaan::vertex(i, H.vertices()[i].nets(), H.vertices()[i].weight()));
            }

            for(auto& n : H.nets()) {
                for(size_t i = 0; i < n.vertices().size(); i++) {
                    n.vertices()[i] = H.local_id(n.vertices()[i]);
                }
            }

            for (size_t i = 0; i < EX; i++) {
                if (H.nets()[i].size() > 0) {
                    netsx.push_back(pmondriaan::net(i, H.nets()[i].vertices(), H.nets()[i].cost()));
                }
            }

            auto HC = pmondriaan::hypergraph(VX, EX, verticesx, netsx, nzx);

            HC(H.local_id(idf)).set_fixpart(fixed_vertex_part);
            HC(H.local_id(idf)).set_fixwdiff(fix_weight_difference);
            HC(H.local_id(idf)).set_fixcdiff(fc);

            std::cout << "lol: " << fc << "\n";


            pmondriaan::remove_free_nets(world, HC, 0);


            
            std::cout << "transform over\n";
            return std::move(HC);
        }


        pmondriaan::remove_free_nets(world, H, 0);

        std::cout << "transform over 2\n";

        return std::move(H);
    }

    // List of fixed vertices of each bisection part.
    auto fixed1 = std::vector<long>();
    auto fixed2 = std::vector<long>();
    std::string ffline;

    // Read fixed file
    while (std::getline(ffin, ffline)) {
        try {
            std::vector<long> pair = std::vector<long>();
            char *pch;
            const char s[2] = " ";
            pch = strtok ((char*)ffline.c_str(),s);
            while (pch != NULL) {
                pair.push_back(std::stol(pch));
                pch = strtok (NULL, " ");
            }
            if(pair[1] == 0) {
                fixed1.push_back(pair[0]);
                H(H.local_id(pair[0])).set_fixpart(0);
            }
            if(pair[1] ==  1){
                fixed2.push_back(pair[0]);
                H(H.local_id(pair[0])).set_fixpart(1);
            }
        } catch (...) { std::cerr << "Error: fixfile\n"; }
    }

    // We need to transform due to fixed vertices in both parts.
    if(fixed1.size() != 0 && fixed2.size() != 0) {
        std::cout << "transform fixed\n";
        auto fix_nets1 = std::vector<long>();
        auto fix_nets2 = std::vector<long>();
        auto fix_nets_overlap = std::vector<long>();
        auto netsf1 = std::vector<long>();
        auto netsf2 = std::vector<long>();
        long w1 = 0;
        long idf1 = fixed1[0];
        long idf2 = fixed2[0];
        long g = 0;
        long maxid = 0;

        long fc = 0;
        long fc1 = 0;
        long fc2 = 0;
        for(auto& n : H.nets()) {
            auto vec = std::vector<long>();
            bool net_has_fix1 = false;
            bool net_has_fix2 = false;
            if(n.id() > maxid) {
                maxid = n.id();
            }
            for(auto v : n.vertices()) {
                if(H(H.local_id(v)).fixpart() == 0) {
                    net_has_fix1 = true;
                } else if(H(H.local_id(v)).fixpart() == 1) {
                    net_has_fix2 = true;
                } else {
                    vec.push_back(v);
                }
            }
            if(net_has_fix1 == true && net_has_fix2 == true) {
                fc += n.cost();
                fix_nets_overlap.push_back(n.id());
            } else if(net_has_fix2 == true) {
                fc2 += n.cost();
                fix_nets2.push_back(n.id());
                netsf2.push_back(n.id());
            } else if(net_has_fix1 == true) {
                fc1 += n.cost();
                fix_nets1.push_back(n.id());
                netsf1.push_back(n.id());
            }

            while(n.vertices().size() > 0) {
                n.pop_back();
            }
            for(auto indi : vec) {
                n.add_vertex(indi);
            }
            
            g++;
        }
        std::cout << "transform fixed1\n";

        for(auto n : fix_nets_overlap) {
            H.remove_net_by_index(H.local_id_net(n));
        }

        for(auto i : fixed1) {
            w1 += H(H.local_id(i)).weight();
            H.remove_free_vertex(i);
            g++;
        }
        std::cout << "transform fixed2\n";

        long w2 = 0;
        for(auto i : fixed2) {
            w2 += H(H.local_id(i)).weight();
            H.remove_free_vertex(i);
            g++;
        }

        std::cout << "transform fixed3 " << H.nets().size() << "\n";

        long mincost = 0;
        long maxcost = 0;
        for(auto& n : H.nets()) {
            if(mincost > n.cost()) {
                mincost = n.cost();
            }
            if(maxcost < n.cost()) {
                maxcost = n.cost();
            }
        }
        std::cout << "transform fixed4 minco " << mincost << "\n";
        std::cout << "transform fixed4 maxco " << maxcost << "\n";

        std::cout << "transform fixed4 " << H.nets().size() << "\n";

        std::cout << "Fix nets: " << fix_nets1.size() + fix_nets2.size() << "\n";

        auto minw = w1;
        if(w1 > w2) {
            minw = w2;
        }

        H.add_vertex(idf1, netsf1, w1 - minw + 1);
        H.add_vertex(idf2, netsf2, w2 - minw + 1);

        H.add_to_nets(H(H.local_id(idf1)));
        H.add_to_nets(H(H.local_id(idf2)));

        auto vertsf = std::vector<long>();
        vertsf.push_back(idf1);
        vertsf.push_back(idf2);

        H.add_net(maxid + 1, vertsf, -1000000);

        auto val = 0;
        for(auto n : H.nets()) {
            val += n.size();
        }

        size_t nzx = val;
        size_t VX = V - fixed1.size() - fixed2.size() + 1;
        size_t EX = H.nets().size();
        auto verticesx = std::vector<pmondriaan::vertex>();
        auto netsx = std::vector<pmondriaan::net>();

        for (size_t i = 0; i < H.vertices().size(); i++) {
            for(size_t j = 0; j < H.vertices()[i].nets().size(); j++) {
                H.vertices()[i].nets()[j] = H.local_id_net(H.vertices()[i].nets()[j]);
            }
            verticesx.push_back(pmondriaan::vertex(i, H.vertices()[i].nets(), H.vertices()[i].weight()));
        }

        for(auto& n : H.nets()) {
            for(size_t i = 0; i < n.vertices().size(); i++) {
                n.vertices()[i] = H.local_id(n.vertices()[i]);
            }
        }

        for (size_t i = 0; i < EX; i++) {
            if (H.nets()[i].size() > 0) {
                netsx.push_back(pmondriaan::net(i, H.nets()[i].vertices(), H.nets()[i].cost()));
            }
        }

        auto HC = pmondriaan::hypergraph(VX, EX, verticesx, netsx, nzx);

        HC(H.local_id(idf1)).set_fixpart(0);
        HC(H.local_id(idf1)).set_fixwdiff(minw - 1);
        HC(H.local_id(idf1)).set_fixcdiff(1000000 + fc);

        std::cout << "lol: " << fc << "\n";


        pmondriaan::remove_free_nets(world, HC, 0);


        
        std::cout << "transform over\n";
        return std::move(HC);
    }


    pmondriaan::remove_free_nets(world, H, 0);

    std::cout << "transform over 2\n";

    return std::move(H);
}

std::optional<pmondriaan::hypergraph> read_hypergraph(std::string file, std::string mode_weight, std::string fixfile) {
    std::cout << fixfile << "\n";
    std::ifstream fs(file);
    if (fs.fail()) {
        std::cerr << "Error: " << std::strerror(errno);
        return std::nullopt;
    }
    std::ifstream ffs(fixfile);
    if (ffs.fail()) {
        std::cerr << "Error: " << std::strerror(errno);
        return std::nullopt;
    }
    return read_hypergraph_istream(fs, ffs, mode_weight);
}

std::optional<pmondriaan::hypergraph>
read_hypergraph(std::string file, bulk::world& world, std::string mode_weight, std::string fixfile) {
    std::cout << fixfile << "\n";
    std::ifstream fs(file);
    if (fs.fail()) {
        std::cerr << "Error: " << std::strerror(errno);
        return std::nullopt;
    }
    std::ifstream ffs(fixfile);
    if (ffs.fail()) {
        std::cerr << "Error: " << std::strerror(errno);
        return std::nullopt;
    }
    return read_hypergraph_istream(fs, world, ffs, mode_weight);
}

} // namespace pmondriaan
