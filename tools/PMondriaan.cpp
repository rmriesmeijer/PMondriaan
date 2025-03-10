#include <iostream>
#include <string>

#include <CLI/CLI.hpp>

#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
using environment = bulk::mpi::environment;
#else
#include <bulk/backends/thread/thread.hpp>
using environment = bulk::thread::environment;
#endif

#include <pmondriaan.hpp>

int main(int argc, char** argv) {

    /* Sequential reading of parameters */
    CLI::App app("PMondriaan settings");
    app.option_defaults()->required();

    app.set_config("--config", "../tools/defaults.toml", "Read a TOML file", true);

    struct cli_settings {
        int p = 2;
        long k = 2;
        double eps = 0.05;
        double eta = 0.10;
        std::string matrix_file;
        std::string hypergraph_weights;
        std::string breaking_mode;
        std::string limit_edge_size;
        std::string simplify_mode;
    };

    cli_settings settings;

    auto options = pmondriaan::options();

    CLI::Option* fopt = app.add_option("-f, --file", settings.matrix_file,
                                       "File including the hypergraph to be "
                                       "partitioned in matrixmarket format");
    fopt->check(CLI::ExistingFile);

    app.add_option("-p, --processors", settings.p,
                   "Number of processors to be used", settings.p);
    app.add_option("-k, --parts", settings.k,
                   "Number of parts to partition H into", settings.k);

    app.add_option("--eps, --epsilon", settings.eps,
                   "Maximum imbalance of the final partitioning", settings.eps);
    app.add_option("--eta", settings.eta,
                   "Maximum imbalance during the parallel computation", settings.eta);

    CLI::Option* wopt =
    app.add_option("--weights", settings.hypergraph_weights,
                   "How the weights of the vertices should be computed");

    CLI::Option* bopt =
    app.add_option("--breaking_mode", settings.breaking_mode,
                   "If and how hyperedges are broken up, none value will skip breaking triples, break_triples_in_initial_partitioning will use triple breaking right before initial partitioning only");

    CLI::Option* lopt =
    app.add_option("--limit_edge_size", settings.limit_edge_size,
                   "If edge sizes are limited in coarsening, this will ignore any vertices beyond the limit in sequential match finding for vertices");

    CLI::Option* sopt =
    app.add_option("--simplify_mode", settings.simplify_mode,
                   "Where to deploy simplification, complete will utilize every method available, parinit will do parallel simplification and initial partitioning but not sequential simplification, sequential only does sequential simplification and initial will only do initial partitioning simplification.");

    std::map<std::string, pmondriaan::bisection> bisection_map{
    {"random", pmondriaan::bisection::random},
    {"multilevel", pmondriaan::bisection::multilevel}};

    app
    .add_option("--bisect", options.bisection_mode, "The bisection mode used")
    ->transform(CLI::CheckedTransformer(bisection_map, CLI::ignore_case));

    std::map<std::string, pmondriaan::sampling> sampling_map{
    {"random", pmondriaan::sampling::random},
    {"label_propagation", pmondriaan::sampling::label_propagation}};

    app
    .add_option("--sampling", options.sampling_mode, "Sampling mode to be used")
    ->transform(CLI::CheckedTransformer(sampling_map, CLI::ignore_case));

    std::map<std::string, pmondriaan::m> metric_map{{"cutnet", pmondriaan::m::cut_net},
                                                    {"lambda_minus_one",
                                                     pmondriaan::m::lambda_minus_one}};
    app.add_option("--metric", options.metric, "Metric to optimized")
    ->transform(CLI::CheckedTransformer(metric_map, CLI::ignore_case));

    wopt->check(CLI::IsMember({"one", "degree"}));

    bopt->check(CLI::IsMember({"none", "break_triples_in_initial_partitioning"}));

    lopt->check(CLI::IsMember({"true", "false"}));

    sopt->check(CLI::IsMember({"complete", "parinit", "sequential", "initial"}));

    app.add_option("--sample_size", options.sample_size, "The sample size used in the coarsening");
    app.add_option("--coarsening_max_edge_size", options.coarsening_max_edge_size, "The max edge size used in the coarsening");
    app.add_option("--max_cluster_size", options.coarsening_max_clustersize,
                   "The maximum clustersize during coarsening");
    app.add_option("--lp_max_iter", options.lp_max_iterations,
                   "The maximum number of iterations in the label propagation");
    app.add_option("--coarsening_nrvertices", options.coarsening_nrvertices,
                   "The number of vertices in the hypergraph when coarsening "
                   "stops");
    app.add_option("--coarsening_max_rounds", options.coarsening_maxrounds,
                   "The maximum number of coarsening rounds");
    app.add_option("--KLFM_max_passes", options.KLFM_max_passes,
                   "The maximum number of passes during the KLFM algorithm");
    app.add_option("--KLFM_max_no_gain_moves", options.KLFM_max_no_gain_moves,
                   "The maximum number of moves that do not give any positive "
                   "gain during the KLFM algorithm");
    app.add_option("--KLFM_par_send_moves", options.KLFM_par_number_send_moves,
                   "The number of moves to find before synchronizing in the "
                   "parallel KLFM algorithm");

    CLI11_PARSE(app, argc, argv);

    environment env;
    //auto timert = bulk::util::timer();
    /* Start parallel part */
    env.spawn(settings.p, [&settings, &options, &app](bulk::world& world) {
        auto s = world.rank();

        //auto time = bulk::util::timer();

        if (s == 0) {
            // write the settings to the defaults file
            auto conf = app.config_to_str();
            std::ofstream out("../tools/settings_run.toml");
            out << conf;
            out.close();
        }
        //world.log(settings.matrix_file);
        auto hypergraph = pmondriaan::read_hypergraph(settings.matrix_file, world,
                                                      settings.hypergraph_weights);

        if (!hypergraph) {
            std::cerr << "Error: failed to load hypergraph\n";
            return;
        }
        auto H = hypergraph.value();

        auto time = bulk::util::timer();
        recursive_bisect(world, H, settings.k, settings.eps, settings.eta, options, settings.breaking_mode, settings.limit_edge_size, settings.simplify_mode);
        auto time_used = time.get();

        auto lb = pmondriaan::load_balance(world, H, settings.k);
        auto cutsize = pmondriaan::cutsize(world, H, options.metric);
        if (!partitioning_to_file(world, H,
                                  "../tools/results/" +
                                  settings.matrix_file.substr(
                                  settings.matrix_file.find_last_of('/') + 1) +
                                  "-k" + std::to_string(settings.k) + "-p" +
                                  std::to_string(settings.p),
                                  settings.k)) {
            std::cerr << "Error: failed to write partitioning to file\n";
            return;
        }

        auto weight_parts = pmondriaan::global_weight_parts(world, H, settings.k);
        if (s == 0) {
            world.log("Partitioned hypergraph with %d vertices into %d parts "
                      "using %d processors",
                      H.global_size(), settings.k, world.active_processors());
            world.log("Load balance of partitioning found: %lf", lb);
            world.log("Cutsize of partitioning found: %d", cutsize);
            world.log("Time used: %lf milliseconds", time_used);

            for (int i = 0; i < settings.k; i++) {
                world.log("Weight part %d: %ld", i, weight_parts[i]);
            }
        }

        world.sync();
	//world.log("Total main time: %lf", time.get_change());
    });
    //std::cout << "Total main time: " << timert.get_change() << "\n";
    return 0;
}
