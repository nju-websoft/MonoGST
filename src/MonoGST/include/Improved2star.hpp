#pragma once
#include "Graph.hpp"
#include "GlobalUtils.hpp"
#include "Log.hpp"
#include <vector>
using namespace std;

template<typename edgetype>
class Improved2star
{
private:
    edgetype INF;
    int n, g;
    Tree <edgetype> answer;
    edgetype now_min_ans;
    vector <vector <edgetype> > dist;
    vector <vector <int> > prev;
    vector <vector<int>> LS;
    vector <vector<int>> rk;
    vector <edgetype> D_1;
    vector <edgetype> Dsum;

    vector <vector<pair<int, edgetype>> > ed;
    vector <int> deg;
    vector <vector <int> > cover_groups;

    double start_dijk;

    double min_full_cover;
    int min_full_cover_v;

    double average_point_ratio_0, average_point_ratio_all;
    vector<int> dijk_rk;
    double dijk_rk_sum;
    int dijk_rk_cnt;
    vector <int> cover_size_cnt;
    vector <double> cover_size_dijk;
    void init_debug_info()
    {
        average_point_ratio_0 = 0, average_point_ratio_all = 0;
        dijk_rk = vector<int> (n + 1, 0);
        dijk_rk_sum = 0;
        dijk_rk_cnt = 0;
        cover_size_cnt = vector <int> (g + 1, 0);
        cover_size_dijk = vector <double> (g + 1, 0);
    }

    int parameter; // 末位为 0 表示启动 val_2

    vector <edgetype> dis;
    vector <int> pre;
    priority_queue <pair<edgetype, int>, vector<pair<edgetype, int>>, greater<pair<edgetype, int>>> pq;
    vector <int> useful_points;
    vector <int> best;
    vector <int> prev_uncover;
    vector <double> sum;
    vector <int> query_cover;
    int res_uncover;
    double next_aver_mx;
    double now_dijk_mx;
    pair<double, int> chs;
    int current_cover;


    /* 预处理 query 集合顺序 和 trivial 情况提前返回 */
    bool init(const Graph<edgetype>& graph, vector<vector<int>>& query)
    {
        answer = Tree<edgetype>(INF);
        n = graph.get_n(), g = query.size();

        if(g == 1) 
        {
            answer = Tree<edgetype>(0);
            return false;
        }

        for(int i = 1; i < g; i++)
            if(query[i].size() < query[0].size())
                swap(query[i], query[0]);
        
        if(query[0].size() == 0)
            return false;

        vector <int> cover_weight(n + 1, 0);
        for(int i = 0; i < g; i++)
        {
            for(auto x : query[i])
                cover_weight[x] ++;
        }
        sort(query[0].begin(), query[0].end(), [&](int i, int j){return cover_weight[i] > cover_weight[j];});
        
        if(cover_weight[query[0][0]] == g)
        {
            answer = Tree<edgetype>(0);
            return false;
        }

        ed = vector <vector<pair<int, edgetype>> > (n + 1);
        deg = vector <int> (n + 1, 0);
        cover_groups = vector <vector <int> > (n + 1);
        for(int i = 0; i < g; i++)
            for(auto x : query[i])
                cover_groups[x].push_back(i);
        
        return true;
    }

    /* 预处理询问组到每个点的最短路 */
    void init_dist(const Graph<edgetype>& graph, const vector<vector<int>>& query)
    {
        dist = vector <vector <edgetype> > (g, vector <edgetype> (n + 1, INF));
        prev = vector <vector <int> > (g, vector <int> (n + 1, -1));
        LS = vector<vector<int>> (n + 1, vector<int>(g));

        for (auto& v : LS) iota(v.begin(), v.end(), 0);
        edgetype max_dis = 0, D = graph.get_min_weight();

        vector <vector <pair <edgetype, int> > > parti_queue;
        int queue_cur = 0;
        auto reset_queue = [&]()
        {
            queue_cur = 0;
        };
        auto enqueue = [&](int u, edgetype d)
        {
            if(d > max_dis) max_dis = d;
            int id = d / D;
            if((int)parti_queue.size() <= id) parti_queue.resize(id + 1);
            parti_queue[id].push_back({d, u});
        };
        auto parti_empty = [&]()
        {
            while(queue_cur < (int)parti_queue.size() && parti_queue[queue_cur].empty()) queue_cur++;
            return queue_cur == (int)parti_queue.size();
        };
        auto dequeue = [&]()
        {
            while(queue_cur < (int)parti_queue.size() && parti_queue[queue_cur].empty()) queue_cur++;
            auto [d, u] = parti_queue[queue_cur].back();
            parti_queue[queue_cur].pop_back();
            return pair<edgetype, int>{d, u};
        };


        for(int ki = 0; ki < g; ki++)
        {
            vector<int> extended(n + 1);

            for(auto x : query[ki])
                dist[ki][x] = 0, enqueue(x, 0);
            while(!parti_empty())
            {
                auto [d, u] = dequeue();
                if(extended[u]) continue;
                extended[u] = 1;
                max_dis = max(max_dis, d);
                for(auto [v, w] : graph.get_adj()[u])
                {
                    if(dist[ki][v] / D > (dist[ki][u] + w) / D) 
                        enqueue(v, w + dist[ki][u]);

                    if(dist[ki][v] > dist[ki][u] + w)
                        dist[ki][v] = dist[ki][u] + w, prev[ki][v] = u;
                }
            }   
            // Log::info("parti_queue.size() : " + to_string(parti_queue.size()));
            reset_queue();
        }


        rk = vector<vector<int>> (n + 1, vector<int>(g));
        Dsum = vector <edgetype> (g + 1, INF);
        D_1 = vector <edgetype> (g + 1, INF);
        for(int c = 1; c <= n; c++)
        {
            sort(LS[c].begin() + 1, LS[c].end(), [&](int i, int j){return dist[i][c] < dist[j][c];});

            edgetype sum=0;
            for(int j = 1; j <= g - 1; j++) {
                rk[c][LS[c][j]] = j;
                sum += dist[LS[c][j]][c];
                Dsum[j] = min(Dsum[j], sum);
                D_1[LS[c][j]] = min(D_1[LS[c][j]], dist[LS[c][j]][c]);
            }
        }
    }
    
    /* 计算简单估值和对应中心 */
    void calc_min_full_cover()
    {
        min_full_cover = numeric_limits<double>::max();;
        min_full_cover_v = 0;
        for(int i = 1; i <= n; i++) 
        {
            edgetype sum = 0;
            for(int j = 0; j <= g - 1; j++) 
                sum += dist[LS[i][j]][i];
            if(sum / (g - 1) < min_full_cover)
            {
                min_full_cover = 1.0*sum / (g - 1);
                min_full_cover_v = i;
            }
        }

        Tree<edgetype> now_answer;
        int c = min_full_cover_v;
        for(int i = 0; i < g; i++)
            add_edge(c, prev[LS[c][i]], dist[LS[c][i]], now_answer);
        answer = now_answer;
        now_min_ans = answer.get_sum_weight();
    }

    /* 新选根时的预处理 */
    void init_now_root()
    {
        dis = vector <edgetype> (n + 1, INF);
        pre = vector <int> (n + 1, -1);
        while(!pq.empty())
            pq.pop();
        useful_points = vector <int> ();
        best = vector <int> (n + 1, g - 1);
        prev_uncover = vector <int> (n + 1, g - 1);
        sum = vector <double> (n + 1, 0);
        query_cover = vector <int> (g, 0);
        query_cover[0] = 1;
        res_uncover = 0;
        now_dijk_mx = INF;
        chs = {numeric_limits<double>::max(), -1};

    }

    /* 更新当前状态的 dijk 最大值 */
    void find_dis_allowed_mx()
    {
        double now_rig_mx = min(chs.first, next_aver_mx);

        if(current_cover == g - 1)
        {
            now_dijk_mx = now_rig_mx - D_1[res_uncover];
            return;
        }
        now_dijk_mx = 0;
        for(int i = 1; i <= g - current_cover; i++)
            now_dijk_mx = max(now_dijk_mx, now_rig_mx*i - Dsum[i]);
    }

    /* 给 now_answer 添加一条边 */
    bool add_edge(int u, vector <int> &par, vector <edgetype> &weight, Tree<edgetype> &now_answer, vector <int> *new_points = nullptr)
    {
        while(par[u] != -1)
        {
            now_answer.add_edge(u, par[u], weight[u] - weight[par[u]]);
            if(now_answer.get_sum_weight() >= answer.get_sum_weight()) 
                return false;
            u = par[u];
            if(new_points) new_points->push_back(u);
        }   
        return true; 
    }

    /* dijk 新确定一个节点 dis(r,u) */
    void add_point(int u)
    {
        useful_points.push_back(u);
        dijk_rk[u] = useful_points.size();

        sum[u] = dis[u];
        prev_uncover[u] = 0;
        best[u] = 0;

        int n_prev_uncover = 0;
        edgetype n_sum = 0;
        for(int j = 1; j <= g - 1; j++)
        {
            int qid = LS[u][j];
            if(!query_cover[qid])
            {
                n_sum = sum[u] + dist[qid][u];
                n_prev_uncover = prev_uncover[u] + 1;

                if(n_prev_uncover == 1 || 1.0*n_sum / n_prev_uncover < 1.0*sum[u] / prev_uncover[u])
                {
                    best[u] = j;
                    sum[u] = n_sum;
                    prev_uncover[u] = n_prev_uncover;
                }
                else
                    break;
            }
        }


        if(1.0*sum[u] / prev_uncover[u] < chs.first)
        {
            chs = {1.0*sum[u] / prev_uncover[u], u};
            find_dis_allowed_mx();
           // cerr << "u = " << u << " new dub = " << now_dijk_mx << endl;
        }
    }
    /* 继续跑一段 dijk */
    void run_dijkstra(const Graph<edgetype>& graph)
    {
        while(!pq.empty())
        {
            if(pq.top().first >= now_dijk_mx)
                return;
            auto [d, u] = pq.top();
            pq.pop();
            if(d > dis[u]) continue;

            add_point(u);

            for(auto [v, w] : graph.get_adj()[u])
            {
                if(dis[v] > dis[u] + w)
                    dis[v] = dis[u] + w, pre[v] = u, pq.push({dis[v], v});
            }
        }
    }

    /* 更新当前状态的 res_uncover */
    void updata_res_uncover()
    {
        if(current_cover == g - 1)
        {
            for(int i = 1; i <= g - 1; i++) if(!query_cover[i])
                res_uncover = i;
        }
    }

    /* 删叶子并尝试更新答案 */
    void Delete_Leave(Tree<edgetype> &now_answer)
    {
        set <int> vex;
        for(auto [uv, w]: now_answer.get_edges())
        {
            auto [u, v] = uv;
            vex.insert(u), vex.insert(v);
            deg[u]++, deg[v]++;
            ed[u].push_back({v, w});
            ed[v].push_back({u, w});
        }
        vector <int> cover_cnt(g, 0);
        for(auto x : vex)
            for(auto y : cover_groups[x])
                cover_cnt[y]++;
        
        for(auto x: cover_cnt)
            if(x == 0)
                return;

        auto ck = [&](int x)
        {
            for(auto y : cover_groups[x])
                if(cover_cnt[y] == 1)
                    return false;
            for(auto y: cover_groups[x])
                cover_cnt[y]--;
            return true;
        };
        queue <int> q;
        for(auto x : vex)
            if(deg[x] == 1 && ck(x))
                q.push(x);
        while(!q.empty())
        {
            int u = q.front();
            q.pop();
            for(auto [v, w] : ed[u])
            {
                deg[v]--;
                now_answer.delete_edge(u, v);
                if(deg[v] == 1 && ck(v))
                    q.push(v);
            }
        }
        
        for(auto x : vex)
        {
            deg[x] = 0;
            ed[x].clear();
        }

        if(now_answer.get_sum_weight() < now_min_ans) 
        {
            now_min_ans = now_answer.get_sum_weight();
            answer = now_answer;
        }
    }

    /* 没有中间点的情况 */
    void no_mid_point()
    {
        vector <int> cover_by_v(g, 0);
        set <int> cover_v;
        int covered_cnt = 0;
        for(auto x : useful_points)
        {
            for(auto y: cover_groups[x])
                if(cover_by_v[y] == 0)
                    cover_by_v[y] = x, cover_v.insert(x), ++covered_cnt;
        }
        if(covered_cnt != g) return;
        Tree<edgetype> now_ans;
        for(auto x : cover_v)
            if(!add_edge(x, pre, dis, now_ans))
                return;
        Delete_Leave(now_ans);
    }

    /* 根据根和第一层点构建最优树 */
    void Build_Tree(vector <int> &transfer_list)
    {
        Tree<edgetype> now_ans;
        vector <int> points = transfer_list;
        for(auto x : transfer_list)
            if(!add_edge(x, pre, dis, now_ans, &points))
                return;
    
        vector <pair <edgetype, int> > mi_dis(g, {INF, -1});
        vector <int> vis(g, 0);

        auto add_v = [&]()
        {
            for(auto v : points)
            {
                for(int i = 1; i < g; i++)
                    mi_dis[i] = min(mi_dis[i], {dist[i][v], v});
            }
            points.clear();
        };

        add_v();

        for(int i = 1; i < g; i++)
        {
            int mij = -1;
            for(int j = 1; j < g; j++)
                if(!vis[j] && (mij == -1 || mi_dis[j] < mi_dis[mij]))
                    mij = j;
            
            vis[mij] = 1;
            if(!add_edge(mi_dis[mij].second, prev[mij], dist[mij], now_ans, &points))
                return;
            add_v();
        }
        
        Delete_Leave(now_ans);
    }

    void Reset_cost(int c)
    {
        sum[c] = 0;
        prev_uncover[c] = 0;
        best[c] = 0;
        for(int j = 1; j <= g - 1; j++)
        {
            int qid = LS[c][j];
            if(!query_cover[qid])
            {
                sum[c] += dist[qid][c];
                prev_uncover[c] += 1;
                best[c] = j;
                break;
            }
        }   
    }

public:
    Improved2star()
    {
        INF = numeric_limits<edgetype>::max();
        answer = Tree<edgetype>(INF);
    }

    void set_parameter(int parameter)
    {
        this->parameter = parameter;
    }
    
    Tree<edgetype> solve(const Graph<edgetype>& graph, vector<vector<int>>& query)
    {
        // Timer::start("init");
        if(!init(graph, query))
            return answer;

        init_dist(graph, query);
        //Log::info("Init dist end");
        calc_min_full_cover();
        //Log::info("Calc min full cover end");

        init_debug_info();

        // Timer::stop("init", LogLevel::LOG_INFO);
        //Log::info("Init end");

        for(auto r: query[0])
        {
            Tree<edgetype> now_answer;
            vector <int> transfer_list;

            init_now_root();
            
            current_cover = 1;
            next_aver_mx = 1.0*now_min_ans / (g - 1);
            
            dis[r] = 0;
            pq.push({0, r}); 
            run_dijkstra(graph);

            int mx_dijk_rk = 0;

            bool fine = 1;
            double now_sum_weight = 0;

            while(current_cover != g) 
            {
                // cerr << "\n\n";
                updata_res_uncover();

                chs = {numeric_limits<double>::max(), -1};
                for(auto u : useful_points)
                {
                    // cerr << u << " " << sum[u] << " " << prev_uncover[u] << " " << best[u] << endl;
                    if(1.0*sum[u] / prev_uncover[u] < chs.first)
                        chs = {1.0*sum[u] / prev_uncover[u], u};
                }
                next_aver_mx = (1.0*now_min_ans - now_sum_weight) / (g - current_cover);
                find_dis_allowed_mx();
                // cerr << "before " << now_dijk_mx << endl;
                run_dijkstra(graph);

                // if(now_sum_weight + chs.first * (g - current_cover) > now_min_ans)
                // {
                //     fine = 0;
                //     break;
                // }

                int c = chs.second;
                transfer_list.push_back(c);
                int las_cover = current_cover;
                mx_dijk_rk = max(mx_dijk_rk, dijk_rk[c]);
                // cerr << "Best candi" << endl;
                // cerr << chs.first << " " << c << endl;
                // cerr << now_dijk_mx << endl;
                // for(auto x : useful_points)
                //     cerr << x << " " << dis[x] << " " << sum[x] << " " << prev_uncover[x] << " " << best[x] << endl;

                // if(!add_edge(c, pre, dis, now_answer))
                // {
                //     fine = 0;
                //     break;
                // }
                int cur_best_c = best[c];
                for(int i = 1; i <= cur_best_c; i++) if(!query_cover[LS[c][i]])
                {
                    current_cover++;
                    int qid = LS[c][i];
                    query_cover[qid] = 1;
                    // if(!add_edge(c, prev[qid], dist[qid], now_answer)) 
                    // {
                    //     fine = 0;
                    //     break;
                    // }
                    for(auto j : useful_points) 
                    {
                        if(rk[j][qid] <= best[j]) 
                        {
                            sum[j] -= dist[qid][j];
                            prev_uncover[j] -= 1;

                            while(best[j] != g - 1) 
                            {
                                best[j]++;
                                int n_qid = LS[j][best[j]];
                                if(query_cover[n_qid]) continue;

                                if(prev_uncover[j] == 0 || (1.0*sum[j] + dist[n_qid][j]) / (prev_uncover[j] + 1) < 1.0*sum[j] / prev_uncover[j])
                                {
                                    sum[j] += dist[n_qid][j];
                                    prev_uncover[j] += 1;
                                }
                                else 
                                {
                                    --best[j];
                                    break;
                                }    
                            }
                        }
                    }
                }

                Reset_cost(c);

                now_sum_weight += 1.0*chs.first*(current_cover - las_cover);
                

                if(fine == 0) break;
                cover_size_cnt[current_cover]++;
                cover_size_dijk[current_cover] += useful_points.size() * 1.0 / n;
            }

            average_point_ratio_all += useful_points.size() * 1.0 / n;
            dijk_rk_sum += mx_dijk_rk * 1.0 / n;
            dijk_rk_cnt++;

            if(!fine) continue;

            no_mid_point();
            // for(auto c: transfer_list)
            //     cerr << c << " ";
            // cerr << endl;
            Build_Tree(transfer_list);
        }
        
        average_point_ratio_all /= query[0].size();
        dijk_rk_sum /= query[0].size();
        // Log::info("average_point_ratio_all: " + to_string(average_point_ratio_all));
        // Log::info("dijk_rk_sum: " + to_string(dijk_rk_sum));
        // Log::info("query[0].size(): " + to_string(query[0].size()));
        return answer;
    }

    double get_average_point_ratio_all()
    {
        return average_point_ratio_all;
    }

    double get_dijk_rk_sum()
    {
        return dijk_rk_sum;
    }
};