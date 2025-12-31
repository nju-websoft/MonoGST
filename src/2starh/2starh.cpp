#include<bits/stdc++.h>
using namespace std;
vector <vector<pair <int,double > > > edg;
vector <vector <vector <int> > > query;
int n;
void Load_Graph(string graphname)
{
    freopen(("data/" + graphname + "/graph.txt").c_str(), "r", stdin);
    int m;
    cin >> n >> m;
    edg.resize(n + 1);
    for(int i = 0; i < m; i++)
    {
        int u, v;double w;
        cin >> u >> v >> w;
        edg[u].push_back({v,w});
        edg[v].push_back({u,w});
    }
}
void Load_Query(string graphname)
{
    freopen(("data/" + graphname + "/query.txt").c_str(), "r", stdin);
    int q;
    cin >> q;
    query.resize(q);
    for(int i = 0; i < q; i++)
    {
        int g;cin>>g;
        query[i].resize(g);
        for (int j = 0; j < g; j++) 
        {
            int k; cin >> k;                 // 该组的大小
            query[i][j].reserve(k);
            for (int t = 0; t < k; t++) 
            {
                int x; cin >> x;             // 该组的每个节点
                query[i][j].push_back(x);
            }
        }
    }
}
double solve(vector <vector <int> > &query)
{
    auto dijk = [&](int s)
    {
        vector <double> dis(n+1, 1e18);
        vector <int> pre(n+1);
        priority_queue <pair <double, int>, vector <pair <double, int> >, greater <pair <double, int> > > pq;
        pq.push({0, s});
        dis[s] = 0;
        while(!pq.empty())
        {
            auto [d, u] = pq.top();
            pq.pop();
            if(dis[u] < d) continue;
            for(auto [v, w] : edg[u])
            {
                if(dis[v] > dis[u] + w)
                {
                    dis[v] = dis[u] + w;
                    pre[v] = u;
                    pq.push({dis[v], v});
                }
            }
        }
        return make_pair(dis,pre);
    };
    
    int g = query.size();
    vector <vector <double> > dis(n+1);
    vector <vector <int> > pre(n+1);
    for(int i = 1; i <= n; i++) tie(dis[i], pre[i]) = dijk(i);
    cerr << "End Dijk" << endl;
    vector <vector <int> > disg(n + 1, vector <int> (g, -1));
    for(int u = 1; u <= n; u++)
    {
        for(int i = 0; i < g; i++)
        {
            for(auto v: query[i])
            {
                if(disg[u][i] == -1 || dis[u][v] < dis[u][disg[u][i]])
                    disg[u][i] = v;
            }
        }
    }
    cerr << "End Disg" << endl; 
    auto Go_Root_r = [&](int r)
    {
        #define cost(x,y) dis[x][disg[x][y]]
        vector <int> uncovered;
        for(int i = 0; i < g; i++) uncovered.push_back(i);
        double Treesum = 0;

        map <pair <int,int>, int> T;
        auto add_edge = [&](int u,int v)
        {
            if(u>v) swap(u,v);
            if(T.find({u,v}) != T.end()) return 0;
            T[{u,v}] = 1;
            return 1;
        };
        auto add_path = [&](int u,int v)
        {
            while(u != v)
            {
                int nxt = pre[v][u];
                if(add_edge(u, nxt)) Treesum += dis[u][nxt];
                u = nxt;
            }
        };

        while(!uncovered.empty())
        {
            double minorm = 1e18;
            int best_v = -1;
            vector <int> best_prefix;
            for(int v = 1; v <= n; v++)
            {
                sort(uncovered.begin(), uncovered.end(), [&](int i, int j)
                {
                    return cost(v,i)*cost(r,j) <= cost(v,j)*cost(r,i);
                });
                double sumv = 0, sumr = 0;
                double minorm_v = 1e18;
                int best_i = 0;
                for(int i = 0; i < uncovered.size(); i++)
                {
                    sumv += cost(v,uncovered[i]);
                    sumr += cost(r,uncovered[i]);
                    double norm = (dis[r][v] + sumv)/sumr;
                    if(norm < minorm_v)
                    {
                        minorm_v = norm;
                        best_i = i;
                    }
                }
                if(minorm_v < minorm)
                {
                    minorm = minorm_v;
                    best_v = v;
                    best_prefix = vector <int> (uncovered.begin(), uncovered.begin() + best_i);
                }
            }
            
            if(minorm > 1e18) return 1e18;

            for(auto i: best_prefix) 
            {
                uncovered.erase(find(uncovered.begin(), uncovered.end(), i));
                //add r -> best_v -> i
                add_path(r, best_v);
                add_path(best_v, disg[best_v][i]);
            }
        }
        return Treesum;
    };

    double result = 1e18;
    for(int r = 1; r <= n; r++)
    {
        cerr << "Root " << r << endl;
        result = min(result, Go_Root_r(r));
        cerr << "End Root " << r << endl;
    }
       
    if(result == 1e18) return -1;
    else return result;
}
int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        cout << "Usage: " << argv[0] << " <graphname>" << endl;
        return 1;
    }
    // ios::sync_with_stdio(0);cin.tie(0);
    string graphname = argv[1];
    Load_Graph(graphname);
    Load_Query(graphname);
    cerr<<"End Load"<<endl;
    freopen(("results/" + graphname + "_2starh_result.txt").c_str(), "w", stdout);
    for(int i = 0; i < query.size(); i++)
    {
        cerr << "Query " << i << ": " << endl;
        auto time_start = chrono::high_resolution_clock::now();
        auto result = solve(query[i]);
        auto time_end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
        cerr << duration/1000.0 << result << endl;
        printf("%.10f %.10f\n", duration/1000.0, result);
    }
    return 0;
}