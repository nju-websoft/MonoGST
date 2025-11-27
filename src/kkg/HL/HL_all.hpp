#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <limits>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include "../Timer.hpp"
#include "../AdjacencyList.hpp"

using std::vector;
using std::pair;
using std::sort;
using std::min;
using std::swap;
using std::cerr;
using std::string;
using std::endl;
using std::istringstream;
using std::ifstream;
using std::ofstream;
namespace fs = std::filesystem;
template<typename NodeType, typename DistanceType>
class HL_all
{
    using AdjList = AdjacencyList<NodeType, DistanceType>;

    public:
        long long LabelSize;
        double cost_time;
        int is_query_path;
        vector <int> Ls;
        vector<vector<NodeType>> Lv;
        vector<vector<DistanceType>> Ld;
        HL_all(AdjList &lis, string filepre, int is_q_p): graph(lis), is_query_path(is_q_p)
        { 
            //Read_Label(filepre);
            Read_Index(filepre);
        }
        HL_all(AdjList &lis, int is_q_p): graph(lis), is_query_path(is_q_p)
        {
            Timer tim;
            tim.start();
            
            LabelSize = 0;
            //std::cerr << "Start Construct HL!" << std::endl;
            n=graph.GetSize();
            inf=std::numeric_limits<DistanceType>::max();

            Lv.assign(n, vector<NodeType>());
            Ld.assign(n, vector<DistanceType>());
            preid.assign(n, vector <NodeType> ());
            qtable.assign(n, inf);

            vector <int> deg(n);
            for(NodeType u=0;u<n;u++)
                for(auto v: graph.GetAllEdges(u))
                    deg[v.first]++;

            vector <NodeType> vid(n);
            for(NodeType i=0;i<n;i++) vid[i]=i;
            sort(vid.begin(),vid.end(),[&](auto x,auto y){return deg[x]>deg[y];});

            vector <NodeType> ranking(n);
            for (int i = 0; i < n; ++i) ranking[vid[i]] = i;

            vector <DistanceType> d(n,inf);
            vector <NodeType> par(n);
            
            vector <NodeType> isr(n);

            int cur = 1;
            
            for(int i=0;i<n;i++)
            {
                if(i == (int)(1.0*n/10.0*cur))
                {
                    std::cerr << cur << "0% has done!" << std::endl;
                    ++cur;
                }

                vector <NodeType> inq;
                NodeType rt=vid[i];
                d[rt]=0;

                for(size_t idx = 0; idx < Lv[rt].size(); ++idx)
                    qtable[Lv[rt][idx]] = Ld[rt][idx];
                int flg = 1;
                int cnt_Lu = 0;

                std::priority_queue <pair <DistanceType, NodeType>, vector <pair <DistanceType, NodeType> >, std::greater <pair <DistanceType, NodeType>> > q;
                q.push({d[rt],rt});
                while(!q.empty())
                {
                    auto [dis, u]=q.top();
                    q.pop();
                    if(d[u]<dis) continue;
                    inq.push_back(u);

                    if(isr[u]) 
                    {
                        flg = 0;
                        continue;
                    }

                    if(flg || Query_Greater(u, d[u]))
                    {
                        Lv[u].push_back(rt);
                        Ld[u].push_back(d[u]);
                        preid[u].push_back(par[u]);
                        ++cnt_Lu;
                        for(auto [v,w]: graph.GetAllEdges(u))
                            if(d[u] + w < d[v])
                            {
                                d[v] = d[u] + w;
                                par[v] = u;
                                q.push({d[v],v});
                            }
                    }
                }
                for(auto u: inq) d[u] = inf;
                for(size_t idx = 0; idx < Lv[rt].size(); ++idx)
                    qtable[Lv[rt][idx]] = inf;
                isr[rt] = 1;
                Ls.push_back(cnt_Lu);
            }
            cost_time = tim.stop();
            for(int i = 0; i < n; i++) LabelSize += Lv[i].size();
            //std::cerr << "HL Contruct Done!" << std::endl;
        }

        bool Query_Greater(NodeType u, DistanceType d)
        {
            for(size_t i = 0; i < Lv[u].size(); ++i) if(qtable[Lv[u][i]] <= d - Ld[u][i]) return false;
            return true;
        }
        
        DistanceType Query(NodeType u, NodeType v)
        {
            DistanceType mindis=inf;
            if(Lv[u].size()>Lv[v].size()) swap(u,v);
            for(size_t i = 0; i < Lv[u].size(); ++i) qtable[Lv[u][i]] = Ld[u][i];
            for(size_t i = 0; i < Lv[v].size(); ++i) if(qtable[Lv[v][i]] != inf) mindis = min(mindis, qtable[Lv[v][i]] + Ld[v][i]);
            for(size_t i = 0; i < Lv[u].size(); ++i) qtable[Lv[u][i]] = inf;
            return mindis;
        }

        vector <NodeType> Query_Path(NodeType u, NodeType v)
        {
            if(Lv[u].size()>Lv[v].size()) swap(u,v);
            DistanceType mindis=inf;
            NodeType h;
            for(size_t i = 0; i < Lv[u].size(); ++i) qtable[Lv[u][i]] = Ld[u][i];
            for(size_t i = 0; i < Lv[v].size(); ++i) if(qtable[Lv[v][i]] != inf && mindis > qtable[Lv[v][i]] + Ld[v][i]) { mindis = qtable[Lv[v][i]] + Ld[v][i]; h = Lv[v][i]; }
            for(size_t i = 0; i < Lv[u].size(); ++i) qtable[Lv[u][i]] = inf;

            vector <NodeType> r1,r2;
            if(mindis == inf) return r1;
            
            auto get_path = [&](int x)
            {
                vector <NodeType> r;
                while(x != h)
                {
                    r.push_back(x);
                    int siz = Lv[x].size();
                    for(int i = 0; i < siz; i++)
                    {
                        if(Lv[x][i] == h)
                        {
                            x = preid[x][i];
                            break;
                        }
                    }
                }
                return r;
            };
            r1 = get_path(u), r2 = get_path(v);
            r1.push_back(h);
            r1.insert(r1.end(),r2.rbegin(),r2.rend());
            return r1;
        }
        
        void Output_Index(string filepre)
        {
            if (!fs::exists(filepre)) 
            {
                std::cout << "In HL_all Output: Directory does not exist. Creating it..." << std::endl;
                fs::create_directories(filepre);
            }
            cerr << "In HL_all Output: start output index\n";
            string filename = filepre + "/index.bin";
            std::ofstream ofs(filename, std::ios::binary);

            ofs.write(reinterpret_cast<const char*>(&n), sizeof (n));
            for(NodeType i = 0; i < n; i++)
            {
                NodeType len = Lv[i].size();
                ofs.write(reinterpret_cast<const char*>(&len), sizeof (len));
                for(NodeType j = 0; j < len; j++)
                {
                    ofs.write(reinterpret_cast<const char*>(&Lv[i][j]), sizeof (Lv[i][j]));
                    ofs.write(reinterpret_cast<const char*>(&Ld[i][j]), sizeof (Ld[i][j]));
                }
                ofs.write(reinterpret_cast<const char*>(preid[i].data()), sizeof(NodeType) * len);
            }
            ofs.close();
        }

        void Read_Index(string filepre)
        {
            //if(filepre.find("DBLP") != string::npos)
            //{
            //    string target_dir = "/home1/yxyang_group/HL/maxcc_graph/NBPCL_dblp_maxcc.txt";
            //    Read_Label(target_dir);
            //    return;
            //}
            
            if (!fs::exists(filepre)) 
            {
                std::cout << "In HL_all Read: Directory does not exist." << std::endl;
                exit(0);
            }
            cerr << "In HL_all: start read index\n";
            string filename = filepre + "/index.bin";
            ifstream ifs(filename, std::ios::binary);
            ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
            inf = std::numeric_limits<DistanceType>::max();
            Lv.resize(n);
            Ld.resize(n);
            preid.resize(n);
            qtable.resize(n, inf);
            for(NodeType i = 0; i < n; i++)
            {
                NodeType len;
                ifs.read(reinterpret_cast<char*>(&len), sizeof(len));
                Lv[i].resize(len);
                Ld[i].resize(len);
                for(NodeType j = 0; j < len; j++)
                {
                    ifs.read(reinterpret_cast<char*>(&Lv[i][j]), sizeof(Lv[i][j]));
                    ifs.read(reinterpret_cast<char*>(&Ld[i][j]), sizeof(Ld[i][j]));
                }
                preid[i].resize(len);
                ifs.read(reinterpret_cast<char*>(preid[i].data()), sizeof(NodeType) * len);
            }
            ifs.close();
            cerr << "In HL_all: end read index\n";
        }
        
        void Output_Label(string filepre)
        {
            ofstream outputFile(filepre);
            outputFile << std::fixed << std::setprecision(15);
            outputFile << n << '\n';
            long long sum = 0;
            int per = 1;
            for(int i = 0; i < n; i++)
            {
                if(i == (int)(1.0*n/10.0*per))
                {
                    cerr << "Finish print "<< per*10 <<"%\n";
                    ++per;
                }
                for(size_t j = 0; j < Lv[i].size(); ++j)
                    outputFile << Lv[i][j] << " " << Ld[i][j] <<  ",";
                outputFile << '\n';
                sum += Lv[i].size();
            }

            cerr << "The row of HBLL is " << n << '\n';
            cerr << "The sum of label is: " << sum << '\n';
            
            outputFile.close();
        }

        void Read_Label(string filepre)
        {
            clear_all();
            //cerr << filepre << "\n";
            ifstream inputFile(filepre);

            if (!inputFile.is_open()) 
                cerr << "Can't open file: " <<  filepre << endl;
            
            cerr << "In HL_all Read: start read label\n";
            Timer tmp_t;
            tmp_t.start();

            inputFile.tie(nullptr);
            string noused;
            inputFile >> n >> noused;
            Lv.assign(n, vector<NodeType>());
            Ld.assign(n, vector<DistanceType>());
            preid.assign(n, vector <NodeType> ());
            qtable.assign(n, inf);
        
            int per = 1;
            for(NodeType i = 0; i < n; i++)
            {
                if(i == (int)(1.0*n/10.0*per))
                {
                    cerr << "Finish read "<< per*10 <<"%\n";
                    ++per;
                }
                NodeType len;
                inputFile >> noused >> len;
                Lv[i].resize(len);
                Ld[i].resize(len);
                preid[i].resize(len);
                for(NodeType j = 0; j < len; j++)
                    inputFile >> Lv[i][j] >> Ld[i][j] >> preid[i][j];
            }
            inputFile.close();
            cerr << "In HL_all Read: end read label\n";
            cerr << "HL_all ReadTime: " + std::to_string(tmp_t.stop()) + " ms.\n";
        }

    private:
        AdjList &graph;
        vector <vector <NodeType> > preid;
        vector<DistanceType> qtable;
        DistanceType inf;
        NodeType n;

        void clear_all()
        {
            n = 0;
            inf = std::numeric_limits<DistanceType>::max();
            Lv.clear();
            Ld.clear();
            qtable.clear();
        }
};