#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <assert.h>

#include "dst/dst.hpp"

using namespace dst;

int n, m, g;
std::vector<std::pair<int, int>> edges;
std::vector<double> weights;
std::vector<int> terms;

bool sort_by_size(std::vector<int> &x, std::vector<int> &y) {
    return x.size() < y.size();
}

int main(int argc, char** argv) {

    if(argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <graphname>" << std::endl;
        return 1;
    }
    std::string graph_file = argv[1];

    /* DATA */

    
    freopen(("data/" + graph_file + "/graph.txt").c_str(), "r", stdin);
    std::cin >> n >> m;
    for (int i = 0; i < m; ++i) {
        int u, v;
        double w;
        std::cin >> u >> v >> w;
        edges.push_back({u, v});
        weights.push_back(w);
        edges.push_back({v, u});
        weights.push_back(w);
        int max = std::max(u, v);
        assert(max <= n);
    }

    // std::vector<std::pair<int, int>> edges{
    //     std::make_pair(0, 1),  std::make_pair(0, 2),  std::make_pair(0, 3),
    //     std::make_pair(0, 4),  std::make_pair(1, 11), std::make_pair(2, 12),
    //     std::make_pair(3, 13), std::make_pair(4, 11), std::make_pair(4, 12),
    //     std::make_pair(4, 13)};
    // std::vector<double> weights{1, 1, 1, 2.5, 2, 2, 2, 1, 1, 1};
    // int root = 0;

    fclose(stdin);

    freopen(("data/" + graph_file + "/query.txt").c_str(), "r", stdin);

    // std::vector<int> terms{11, 12, 13};
    std::vector<std::vector<int>> groups;
    int T; std::cin >> T;
    for (int cas = 1; cas <= T; ++cas) {
        std::cerr << cas <<" start"<< std::endl;
		std::cin >> g;
		std::cerr << cas <<" g=" << g << std::endl;
        groups.clear();
        groups.resize(g + 1);
        for (int i = 1; i <= g; ++i) {
            int size;
            std::cin >> size;
            for (int j = 1; j <= size; ++j) {
                int u;
                std::cin >> u;
                groups[i].push_back(u);
            }
        }
        std::sort(groups.begin(), groups.end(), sort_by_size);
        for (int i = 2; i <= g; ++i) {
            int new_node = i + n - 1;
            for (auto u : groups[i]) {
                edges.push_back({u, new_node});
                weights.push_back(0);
            } 
            terms.push_back(new_node);
        }

        /* START RUNNING */
        std::string method = "fast_level2";
        double alpha = 1;
        int nthresholds = 50;
        double ans = 1e18;
        std::clock_t c_start = std::clock();
        for (auto root : groups[1]) {
            // one extra run of dijkstra in constructor to remove unreachable vertices
            DST dt = DST(edges, weights, root, terms);
            std::shared_ptr<PartialTreeManager> partree = nullptr;
            if (method.compare("level2") == 0) {
                partree = dt.level2_alg();
            } else if (method.compare("fast_level2") == 0) {
                partree = dt.fast_level2_alg();
            } else if (method.compare("level3") == 0) {
                partree = dt.level3_alg();
            } else if (method.compare("fast_level3") == 0) {
                partree = dt.fast_level3_alg(alpha, nthresholds);
            } else {
                std::cerr << "unknown method: " << method << std::endl;
                return 1;
            }
            std::shared_ptr<Tree> tree = partree->to_tree();
            ans = std::min(ans, tree->cost_trimmed());
        }
        std::clock_t c_end = std::clock();
        freopen(("results/" + graph_file + "_DST_result.txt").c_str(), "w", stdout);
        double time_elapsed_s = (c_end - c_start) / CLOCKS_PER_SEC;
        std::cout << time_elapsed_s << ' ' << ans << std::endl;
        fclose(stdout);

        for (int i = 2; i <= g; ++i) {
            for (auto _ : groups[i]) {
                edges.pop_back();
                weights.pop_back();
            } 
            terms.pop_back();
        }
    }    

    fclose(stdin);




    

    return 0;
}