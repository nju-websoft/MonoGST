#include "../HL/HL_all.hpp"
#include "../KKG/Key_KG_plus.hpp"
#include "../AdjacencyList.hpp"
#include "../Utilities.hpp"
#include "../DistanceOracle.hpp"
#include "../Load_Graph.hpp"
#include "../Timer.hpp"
#include "../Prepare_Test.hpp"
#include "../Load_Query.hpp"

#include <type_traits>
#include <format>
#include <iostream>
#include <sstream>
#include <iomanip>
namespace fs = std::filesystem;
int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <graphname>" << std::endl;
        return 1;
    }
    Parameter para;
    para.Load_Parameters();

    string graphname = argv[1];
    string graph_type = "1";
    
    FileManager fop("./results");
    fop.clearFile(graphname + "_KKG_result.txt");

    auto GoKKG = [&]<typename EdgeType>()
    {
        LoadGraph<int, EdgeType> g("./data/" + graphname, graph_type == "0" ? 0 : 1, 1);
        DistanceOracle<int, EdgeType, HL_all> DO(g.adjList, "KeyKG_Index/" + graphname + "/", 1);
        KeyKG_plus<int, EdgeType, HL_all> kkg(g.adjList, DO);
    
        auto read_query_file_at = [&](const fs::path& query_file_path)
        {
            vector<vector<vector<int>>> all_queries;
            ifstream fin(query_file_path.string());
            int q; fin >> q;
            for (int i = 0; i < q; ++i)
            {
                int g; fin >> g;
                vector<vector<int>> query;
                for (int j = 0; j < g; ++j)
                {
                    int s,v; fin >> s;
                    set <int> group;
                    for (int k = 0; k < s; ++k)
                    {
                        fin >> v;
                        group.insert(v - 1);
                    }
                    query.push_back(vector<int>(group.begin(), group.end()));
                }
                all_queries.push_back(query);
            }
            return all_queries;
        };
    
        auto queries = read_query_file_at("./data/" + graphname + "/query.txt");
        cerr << queries.size() << '\n';
        for(auto qry: queries)
        {
            auto [Tw, Tv, Te] = kkg.Query(qry);
            cerr << to_string(kkg.cost_time/1000.0) + " " + to_string(Tw) << '\n';
            fop.writeToFile(graphname + "_KKG_result.txt",to_string(kkg.cost_time/1000.0) + " " + to_string(Tw));
            if(graph_type == "0")
            {
                string vex;
                for(auto x: Tv) vex += to_string(x) + " ";
                fop.writeToFile(graphname + "_KKG_result.txt", vex);
            }
        }
    };

    if(graph_type == "0") 
        GoKKG.template operator()<int>();
    else
        GoKKG.template operator()<long double>();


    
    return 0;
}