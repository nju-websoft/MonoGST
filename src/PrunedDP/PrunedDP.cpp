/*
pinkyhead 2025/9/15
*/
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <queue>
#include <vector>
#include <assert.h>
#include <memory.h>
#include <cmath>
#include <limits>
#include <chrono>
#define all(x) x.begin(), x.end()
typedef long long ll;
// typedef double db;
typedef long double db;
typedef std::pair<int, int> pii;
typedef std::pair<ll, int> pli;
typedef std::pair<ll, ll> pll;
typedef std::pair<int, db> pid;

struct node {
    int v;
    int X;
    db cost;
    db lower_bound;
    node(int v, int X, db cost, db lower_bound) {
        this->v = v;
        this->X = X;
        this->cost = cost;
        this->lower_bound = lower_bound;
    }
    bool operator<(const node& other) const {
        if (lower_bound == other.lower_bound) return cost > other.cost;
        return lower_bound > other.lower_bound;
    }
};

const db eps = 1e-10;
const db inf = std::numeric_limits<double>::max() / 3;
int n, m, group_count, vertex_offset;
int new_n, mask; // new_n = n + group_count, mask = 1 << group_count
std::vector<std::vector<int>> groups;
std::vector<std::vector<int>> color;
std::vector<std::vector<std::pair<int, db>>> e;
std::vector<std::vector<std::vector<db>>> W_pair;
std::vector<std::vector<db>> W_single;
std::vector<std::vector<std::pair<int, db>>> D;
std::priority_queue<node> Q;
std::vector<std::vector<db>> dist; // dist[u][v]  u : virtual, v : any vertex
std::vector<std::vector<db>> predict;

db best;
std::string graph_file;

namespace ALLPaths {
    struct path {
        int u, v;
        int X;
        db cost;
        path(int u, int v, int X, db cost) {
            this->u = u;
            this->v = v;
            this->X = X;
            this->cost = cost;
        }
        bool operator<(const path& other) const {
            return cost > other.cost;
        }
    };

    std::vector<std::vector<std::vector<int>>> D; 
    std::priority_queue<path> Q;

    void init() {
        D.clear();
        D.resize(group_count + 1);
        for (int i = 1; i <= group_count; ++i) {
            D[i].resize(group_count + 1);
            for (int j = 1; j <= group_count; ++j) {
                D[i][j].resize(mask + 1, 0);
            }
        }
    }
    void cal_W() {
        init();
        for (int i = 1; i <= group_count; ++i) {
            for (int j = 1; j <= group_count; ++j) {
                Q.push(path(i, j, (1 << (i - 1)) | (1 << (j - 1)), dist[i][j + n]));
                Q.push(path(i, j, 1 << (i - 1), dist[i][j + n]));
                Q.push(path(i, j, 1 << (j - 1), dist[i][j + n]));
                Q.push(path(i, j, 0, dist[i][j + n]));
            }
        }
        while (!Q.empty()) {
            auto it = Q.top();
            Q.pop();
            if (D[it.u][it.v][it.X]) continue;
            D[it.u][it.v][it.X] = 1;
            W_single[it.u][it.X] = std::min(W_single[it.u][it.X], W_pair[it.u][it.v][it.X]);
            for (int i = 1; i <= group_count; ++i) {
                int new_X = it.X | (1 << (i - 1));
                db new_cost = it.cost + dist[it.v][i + n];
                if (new_cost < W_pair[it.u][i][new_X]) {
                    W_pair[it.u][i][new_X] = new_cost;
                    Q.push(path(it.u, i, new_X, new_cost));
                }
                if (new_cost < W_pair[it.u][i][it.X]) {
                    W_pair[it.u][i][it.X] = new_cost;
                    Q.push(path(it.u, i, it.X, new_cost));
                }
            }
        }
        for (int i = 1; i <= n; ++i) {
            predict[i][0] = 0;
            for (int j = 1; j <= mask; ++j) {
                db min_1 = inf, max_2 = 0, max_3 = 0;
                for (int u = 1; u <= group_count; ++u) {
                    db min_2 = inf;
                    if (!(j & (1 << (u - 1)))) continue;
                    for (int v = 1; v <= group_count; ++v) {
                        if (!(j & (1 << (v - 1)))) continue;
                        min_1 = std::min(min_1, dist[u][i] + W_pair[u][v][j] + dist[v][i]);
                        min_2 = std::min(min_2, dist[v][i]);
                    }
                    max_2 = std::max(max_2, dist[u][i] + W_single[u][j] + min_2);
                }
                for (int k = 1; k <= group_count; ++k) {
                    if (!(j & (1 << (k - 1)))) continue;
                    max_3 = std::max(max_3, dist[k][i]);
                }
                predict[i][j] = std::max(std::max(min_1 / 2, max_2 / 2), max_3);
            }
        }
    }
}

void cal_dist() {
    std::priority_queue<std::pair<db, int>, std::vector<std::pair<db, int>>, std::greater<std::pair<db, int>>> q;
    for (int i = 1; i <= group_count; ++i) {
        std::vector<bool> used(new_n + 1, 0);
        dist[i][i + n] = 0;
        q.push({0, i + n});
        while (!q.empty()) {
            auto [val, u] = q.top(); q.pop();
            if (used[u]) continue;
            used[u] = 1;
            for (auto [v, w] : e[u]) {
                if (dist[i][v] > dist[i][u] + w) {
                    dist[i][v] = dist[i][u] + w;
                    if (v <= n) { // not virtual vertex
                        q.push({dist[i][v], v});
                    }
                }
            }
        }
    }
    e.resize(n + 1);
    for (int i = 1; i <= n; ++i) {
        while (!e[i].empty() && e[i].back().first > n) e[i].pop_back();
    }
} 

db lb(int v, int X, db cost) {
    int inv_X = X ^ mask;
    return predict[v][inv_X] + cost;
}

void update(int v, int X, db cost, db lower_bound) {
    if (D[v][X].first) {
        //assert(cost >= D[v][X].second);
        return;
    }
    lower_bound = std::max(lb(v, X, cost), lower_bound);
    if (lower_bound >= best) return;
    if (X == mask) best = std::min(best, cost);
    Q.push(node(v, X, cost, lower_bound));
}

db pruned_dp_plus_plus() {
    cal_dist();
    ALLPaths::cal_W();
    for (int i = 1; i <= n; ++i) {
        int X = 0;
        for (auto it : color[i]) {
            X = 1 << (it - 1);
            Q.push(node(i, X, 0, lb(i, X, 0)));
        }
        // Q.push(node(i, 0, 0, lb(i, 0, 0)));
    }
    while (!Q.empty()) {
        auto cur = Q.top(); Q.pop();
        if (D[cur.v][cur.X].first) continue;
        D[cur.v][cur.X] = {1, cur.cost};
        if (cur.cost >= best) continue;
        if (cur.X == mask)  {
            best = cur.cost;
            continue;
        }

        int inv_X = mask ^ cur.X;
        db expect_cost = cur.cost;
        for (int i = 1; i <= group_count; ++i) {
            if ((1 << (i - 1)) & inv_X) {
                expect_cost += dist[i][cur.v];
                expect_cost = std::min(expect_cost, inf);
            }
        }
        best = std::min(best, expect_cost);
        if (D[cur.v][inv_X].first) {
            update(cur.v, mask, cur.cost + D[cur.v][inv_X].second, cur.lower_bound);
        }
        if (cur.cost <= best / 2 + eps) {
            for (auto [u, val] : e[cur.v]) {
                assert(u <= n);
                db new_cost = cur.cost + val;
                update(u, cur.X, new_cost, cur.lower_bound);
            }
            for (int i = inv_X; i >= 1; i = (i - 1) & inv_X) {
                if (!D[cur.v][i].first) continue;
                db new_cost = cur.cost + D[cur.v][i].second;
                if (new_cost <= 2.0 / 3 * best + eps) {
                    update(cur.v, i | cur.X, new_cost, cur.lower_bound);
                }
            }
        }
    }
    return best;
}

void clear_all_stuff() {
    new_n = n + group_count;
    groups.clear();
    groups.resize(group_count + 1);
    color.clear();
    color.resize(n + 1);
    while (!Q.empty()) Q.pop();
    for (int i = 1; i <= n; ++i) {
        while (!e[i].empty()&&e[i].back().first > n) e[i].pop_back();
    }
    mask = (1 << group_count) - 1;
    W_pair.clear();
    W_pair.resize(group_count + 1);
    for (int i = 1; i <= group_count; ++i) {
        W_pair[i].resize(group_count + 1);
        for (int j = 1; j <= group_count; ++j) {
            W_pair[i][j].resize(mask + 1, inf);
        }
    }
    W_single.clear();
    W_single.resize(group_count + 1);
    for (int i = 1; i <= group_count; ++i) {
        W_single[i].resize(mask + 1, inf);
    }
    D.clear();
    D.resize(n + 1);
    for (int i = 1; i <= n; ++i) {
        D[i].resize(mask + 1, {0, inf});
    }

    best = inf;
    // W, D, e (virtual nodes), Q, best
    dist.clear();
    dist.resize(group_count + 1);
    for (int i = 1; i <= group_count; ++i) {
        dist[i].resize(new_n + 1, inf);
    }
    predict.clear();
    predict.resize(n + 1);
    for (int i = 1; i <= n; ++i) {
        predict[i].resize(mask + 1, 0);
    }
}

int main(int argc, char* argv[]) {
    if(argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <graphname>" << std::endl;
        return 1;
    }
    std::cerr << "Start reading graph" << std::endl;
    graph_file = argv[1];
    freopen(("data/" + graph_file + "/graph.txt").c_str(), "r", stdin);
    std::cin >> n >> m;
    e.resize(n + 1);
    for (int i = 1; i <= m; ++i) {
        int u, v;
        db w;
        std::cin >> u >> v >> w;
        e[u].push_back({v, w});
        e[v].push_back({u, w});
        if (!u || !v) vertex_offset = 1;
    }
    if (vertex_offset) {
        for (int i = n; i >= 0; --i) {
            e[i] = e[i - 1];
            for (auto &[v, _] : e[i]) {
                v += vertex_offset;
            }
        }
        e[0].clear();
    }
    fclose(stdin);
    std::cin.clear();
    std::cerr << "Finish reading graph" << std::endl;
    freopen(("data/" + graph_file + "/query.txt").c_str(), "r", stdin);
    freopen(("results/" + graph_file + "_PrunedDP_result.txt").c_str(), "w", stdout);
    int T;
    std::cin >> T;
    for (int i = 1; i <= T; ++i) {
        std::cin >> group_count;
        clear_all_stuff();
        for (int j = 1; j <= group_count; ++j) {
            int current_size;
            std::cin >> current_size;
            while (current_size--) {
                int u;
                std::cin >> u;
                groups[j].push_back(u);
                color[u].push_back(j);
            }
        }
        auto st = std::chrono::high_resolution_clock::now();
        e.resize(new_n + 1);
        for (int i = 1; i <= group_count; ++i) {
            int virtual_vertex = i + n;
            for (auto u : groups[i]) {
                e[virtual_vertex].push_back({u, 0});
                e[u].push_back({virtual_vertex, 0});
            }
        }
        db ans = pruned_dp_plus_plus();
        ans = (ans >= inf ? -1 : ans);
        double time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - st).count()/1000000.0;
        std::cout << time << ' ' << ans << '\n';
        std::cerr << time << ' ' << ans << '\n';
    } 
    

    return 0;
}
