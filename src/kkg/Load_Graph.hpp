#ifndef LOADGRAPH_H
#define LOADGRAPH_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <queue>
#include "AdjacencyList.hpp"
using std::vector;
using std::queue;
using std::max;
template<typename NodeType, typename DistanceType>
class LoadGraph
{
    public:
        int N;
        AdjacencyList<NodeType, DistanceType> adjList;
        LoadGraph(std::string filepath, int isweight, int offset)
        {
            load_ini_graph(filepath, isweight, offset);
        }

        LoadGraph(const std::string& filepath, int type, vector <string> &str_info, vector <int> &i_info)
        {
            if(type == -1)
            {
                N = 0;
                return;
            }    
            else if(type == 0)
                load_ini_graph(filepath, i_info[0], i_info[1]);
            else if(type == 1)
                load_DIMACS_graph(filepath, str_info[0], i_info[0]);
            //test_trick();
        }

        void test_trick()
        {
            vector <NodeType> deg(N);
            cerr << N << '\n';
            queue <NodeType> q;
            for(NodeType i = 0; i < N; i++) 
            {
                deg[i] = adjList.GetDeg(i);
                if(deg[i] <= 1) q.push(i);
            }    
            cerr << "Qsize: " << q.size() << '\n';
            NodeType cnt = 0;
            while(!q.empty())
            {
                ++cnt;
                NodeType now = q.front();
                q.pop();
                for(auto [v,_]: adjList.GetAllEdges(now))
                {
                    --deg[v];
                    if(deg[v] == 1)
                        q.push(v);
                }
            }
            cerr << "Graph nodes cnt is " << N << '\n';
            cerr << "Leave cnt is " << cnt << '\n';
        }   


        void load_ini_graph(std::string filepath, int isweight, int offset)
        {
            std::string name[2] = {"/graph.txt","/graph.txt"};
            std::string filename = filepath + name[isweight];
            
            std::ifstream inputFile(filename);

            if (!inputFile.is_open()) 
                std::cerr << "Can't open file: " << filename << std::endl;

            std::string line;
            getline(inputFile, line);

            N = extractIntegers(line)[0];
            adjList = AdjacencyList<NodeType, DistanceType>(N);

            int miv = 1e9, mxv = 0;

            while (getline(inputFile, line)) 
            {
                if(isweight == 0)
                {
                    auto info = extractIntegers(line);
                    if(info[0] == info[1]) continue;
                    adjList.AddUndirectedEdge(info[0] - offset, info[1] - offset, 1);
                    miv = min(miv, info[0]);    
                    miv = min(miv, info[1]);
                    mxv = max(mxv, info[0]);
                    mxv = max(mxv, info[1]);
                }
                else
                {
                    std::istringstream iss(line);
                    int u,v; double w;
                    iss>>u>>v>>w;
                    if(u == v) continue;
                    u -= offset, v -= offset;
                    adjList.AddUndirectedEdge(u, v, w);
                    miv = min(miv, u);
                    miv = min(miv, v);
                    mxv = max(mxv, u);
                    mxv = max(mxv, v);
                }
            }
            //cerr << "Min id is " << miv << " Max id is " << mxv << '\n';
            inputFile.close();
        }

        

        /**
         * @param filepath  图文件所在目录（不以 '\\' 或 '/' 结尾）
         * @param filename  像 "USA"、"NY"、"CAL" 这样的区域标识
         * @param type      0 => 读取 "USA-road-d.<filename>.gr"
         *                  1 => 读取 "USA-road-t.<filename>.gr"
         */
        void load_DIMACS_graph(const std::string& filepath, const std::string& filename, int type)
        {
            const char prefix = (type == 0 ? 'd' : 't');
            std::string target = "USA-road-";
            target += prefix;
            target += ".";
            target += filename;
            target += ".gr";

            std::string fullpath = filepath + "/" + target;

            std::ifstream in(fullpath);
            if (!in.is_open())
                throw std::runtime_error("无法打开图文件: " + fullpath);

            std::string line;
            bool initialized = false;
            while (std::getline(in, line)) 
            {
                if (line.empty() || line[0] == 'c') 
                    continue;
                else if (line[0] == 'p') 
                {
                    // p sp <nodeCount> <arcCount>
                    std::istringstream iss(line);
                    std::string p, sp;
                    int nodeCount, arcCount;
                    iss >> p >> sp >> nodeCount >> arcCount;

                    N = nodeCount;
                    adjList = AdjacencyList<NodeType, DistanceType>(N);
                    initialized = true;
                }
                else if (line[0] == 'a') 
                {
                    if (!initialized) 
                        throw std::runtime_error("在读取任何顶点数之前遇到边定义");
                    // a <u> <v> <w>
                    std::istringstream iss(line);
                    char a;
                    NodeType u, v;
                    DistanceType w;
                    iss >> a >> u >> v >> w;
                    --u, --v;
                    adjList.AddUndirectedEdge(u, v, w);
                }
            }

            if (!initialized)
                throw std::runtime_error("未在文件中找到 'p sp ...' 节点定义");
        }

        void Output_Graph(const std::string& file)
        {
            std::ofstream outFile(file, std::ios::app);
            cerr << file << '\n';
            if (!outFile) 
            {
                std::cerr << "Failed to open the file." << std::endl;
                return;
            }
            int sz = adjList.GetSize();
            outFile << sz << ' ' << adjList.GetEdgeCount() << '\n';
            for(int u = 0; u < sz; u++)
                for(auto [v, w]: adjList.GetAllEdges(u))
                    if(u < v)
                        outFile << u + 1 << ' ' << v + 1 << ' ' << w << '\n';
        }
    
    private:
        std::vector<int> extractIntegers(const std::string& input) 
        {
            std::vector<int> integers;
            std::istringstream iss(input);
            std::string token;

            while (iss >> token) 
            {
                try 
                {
                    int number = stoi(token);
                    integers.push_back(number);
                } catch (const std::invalid_argument& e) 
                {
                    std::cerr << "not a number" << std::endl;
                }
            }

            return integers;
        }
};

#endif