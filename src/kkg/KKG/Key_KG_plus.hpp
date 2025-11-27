#include "../DistanceOracle.hpp"
#include "../AdjacencyList.hpp"
#include "../Timer.hpp"
#include <algorithm>
#include <map>
#include <set>
using std::vector;
using std::cerr;
using std::pair;
using std::tuple;
using std::set;
using std::map;
using std::make_pair;
using std::make_tuple;
template<typename NodeType, typename DistanceType, template<typename, typename> class OracleType>
class KeyKG_plus
{
    public:

    using AdjList = AdjacencyList<NodeType, DistanceType>;
    AdjList &graph;
    DistanceType inf;
    double cost_time;
    DistanceOracle <NodeType, DistanceType, OracleType> &oracle;
    vector < vector <pair <NodeType, DistanceType> > > M;
    KeyKG_plus(AdjList &lis, DistanceOracle<NodeType, DistanceType, OracleType> &DO): graph(lis), oracle(DO) 
    {
        inf = std::numeric_limits<DistanceType>::max(); 
        NodeType N = graph.GetSize();
        M = vector < vector <pair <NodeType, DistanceType> > > (11, vector <pair <NodeType, DistanceType> > (N, {N + 1, inf}));
    }

    tuple <DistanceType, vector <NodeType> ,vector <pair <NodeType, NodeType> > > Query(vector <vector <NodeType> > &K)
    {
        // cerr << "Query: ";
        // for(auto Ki: K)
        // {
        //     for(auto x: Ki)
        //         cerr << x << ' ';
        //     cerr << '\n';
        // }
        // cerr << "\n";
        
        Timer tim;
        tim.start();

        int g = K.size();
        DistanceType miw = inf;
        NodeType N = graph.GetSize();
        vector <vector <NodeType> > revMj(g + 1, vector <NodeType> ());
        while(M.size() <= g + 1) M.push_back(vector <pair <NodeType, DistanceType> > (N, {N + 1, inf}));

        for(int i = 1; i < g; i++)
        {
            for(auto v: K[i])
            {
                for(size_t idx=0; idx<oracle.oracle.Lv[v].size(); ++idx)
                {
                    auto r = oracle.oracle.Lv[v][idx];
                    auto d = oracle.oracle.Ld[v][idx];
                    if(M[i][r].second > d)
                        M[i][r] = {v, d}, revMj[i].push_back(r);
                }
            }
        }

        vector <NodeType> Ux;
        for(int i = 0; i < K[0].size();i++)
        {
            int v0 = K[0][i];
            vector <NodeType> v(g);
            v[0] = v0;
            DistanceType sw = 0;
            for(int j = 1; j < g; j++)
            {
                DistanceType w = inf, nw;
                // for(auto x: K[j])
                // {
                //     nw = oracle.Query_Distance(v0, x);
                //     if(nw < w)
                //     {
                //         w = nw;
                //         v[j] = x;
                //     }
                // }
                for(size_t idx=0; idx<oracle.oracle.Lv[v0].size(); ++idx)
                {
                    auto h = oracle.oracle.Lv[v0][idx];
                    auto dish = oracle.oracle.Ld[v0][idx];
                    auto [y, nw] = M[j][h];
                    if(nw != inf)
                        nw += dish;
                    if(nw < w)
                    {
                        w = nw;
                        v[j] = y;
                    }
                }

                if(w == inf)
                {
                    sw = inf;
                    break;
                }
                sw += w;
            }
            if(sw < miw) 
            {
                miw = sw;
                Ux = v;
            }
        }

        // cerr << "In Ux: ";
        // for(auto u: Ux)
        //     cerr << u << ' ';
        // cerr << '\n';

        vector <NodeType> vex;
        vector <pair <NodeType, NodeType> > edg;
        if(miw == inf) 
        {
            for(int i = 1; i < g; i++)
            {
                for(auto r: revMj[i])
                    M[i][r] = {N, inf};
                revMj[i].clear();
            }
            return make_tuple(-1, vex, edg);
        }
        
        DistanceType mit = inf;
        vector <NodeType> revM;
        for(auto u: Ux)
        {
            auto ins_M0 = [&](int v)
            {
                for(size_t idx=0; idx<oracle.oracle.Lv[v].size(); ++idx)
                {
                    auto r = oracle.oracle.Lv[v][idx];
                    auto d = oracle.oracle.Ld[v][idx];
                    if(M[0][r].second > d)
                        M[0][r] = {v, d}, revM.push_back(r);
                }
            };

            DistanceType sumt = 0;
            set <NodeType> Vt, Ur;
            set <pair <NodeType, NodeType>> edge;
            Vt.insert(u);
            ins_M0(u);
            for(auto x: Ux) if(u != x) Ur.insert(x);
            
            while(!Ur.empty())
            {
                NodeType smi, tmi;
                DistanceType wmi = inf;
                for(auto t: Ur)
                {
                    for(size_t idx=0; idx<oracle.oracle.Lv[t].size(); ++idx)
                    {
                        auto h = oracle.oracle.Lv[t][idx];
                        auto dish = oracle.oracle.Ld[t][idx];
                        auto [s, wnow] = M[0][h];
                        if(wnow != inf)
                            wnow += dish;
                        if(wnow < wmi)
                        {
                            wmi = wnow;
                            smi = s, tmi = t;
                            //cerr << "Dist of " << s << " and " << t << " is " << wnow << '\n';
                        }
                    }
                }
                // for(auto s: Ur)
                //     for(auto t: Vt)
                //     {
                //         wnow = oracle.Query_Distance(s, t);
                //         //cerr << "Dist of " << s << " and " << t << " is " << wnow << '\n';
                //         if(wnow < wmi)
                //         {
                //             wmi = wnow;
                //             smi = s, tmi = t;
                //         }
                //     }
                if(wmi == inf) continue;

                auto p = oracle.Query_Path(smi, tmi);

                //cerr << "[s, t] = [" << smi << ", " << tmi << "]\nPath: ";
                // for(auto x: p) cerr << x << ' ';
                //cerr << '\n';

                DistanceType nowt = 0;

                for(int i = 1; i < p.size(); i++)
                {
                    NodeType u = p[i-1], v = p[i];
                    if(u > v) swap(u, v);
                    if(edge.contains(make_pair(u, v))) continue;
                    edge.insert(make_pair(u, v));
                    nowt += graph.GetEdgeuv(u, v);
                }
                sumt += nowt;
                //cerr << "nowt: " << nowt << '\n';
                for(auto x: p)
                {
                    // cerr << x << ' ';
                    auto it = Ur.find(x);
                    if(it != Ur.end())
                        Ur.erase(it);
                    if(!Vt.contains(x))
                    {
                        ins_M0(x);
                        Vt.insert(x);
                    }
                }

                //cerr << "End Round\n";
            }

            if(sumt < mit)
            {
                mit = sumt;
                vex.clear(), edg.clear();
                for(auto u: Vt) vex.push_back(u);
                for(auto uv: edge)
                    edg.push_back(uv);
            }

            for(auto r: revM) M[0][r] = {N + 1, inf};
            revM.clear();
        }


        for(int i = 1; i < g; i++)
        {
            for(auto r: revMj[i])
                M[i][r] = {N, inf};
            revMj[i].clear();
        }
        cost_time = tim.stop();

        auto ck_ans = [&]()
        {
            if(mit == inf) 
            {
                if(vex.size() == 0) return 0;
                else return 4;
            }
            NodeType ns = vex.size();
            map <NodeType, NodeType> v2id;
            NodeType cnt = 0;
            for(auto x: vex) 
            {
                if(v2id.contains(x)) return 1;//contains repeat nodes
                v2id[x] = ++cnt;
            }    
            vector <NodeType> f(ns + 1);
            for(NodeType i = 1; i <= ns; i++) f[i] = i;
            auto Find = [&](auto &&Find, NodeType x) -> NodeType
            {
                return f[x] = f[x] == x ? x: Find(Find, f[x]);
            };
            for(auto [u, v]: edg)
            {
                auto fu = Find(Find, v2id[u]), fv = Find(Find, v2id[v]);
                f[fv] = fu;
            }
            for(NodeType i = 2; i <= ns; i++) if(Find(Find, i) != Find(Find, 1)) return 2;//Not a tree

            for(auto Ki: K)
            {
                int flg = 0;
                for(auto x: Ki)
                    if(v2id.contains(x))
                        flg = 1;
                if(flg == 0) return 3;//Exist Ki \cap T = \empty
            }
            return 0;//ok
        };

        int cka = ck_ans();
        if(cka != 0)
        {
            cerr << "KeyKG_plus error code is " << ck_ans() << '\n';
        }
        // else
        //     cerr << "Check KeyKG_plus legality ok\n";

        return make_tuple(mit, vex, edg);
    }
};