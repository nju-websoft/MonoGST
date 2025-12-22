#include "../include/Graph.hpp"
#include "../include/Log.hpp"
#include "../include/GlobalUtils.hpp"
#include "../include/Improved2star.hpp"
#include "../include/Improved2star_Abl.hpp"
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <graphname> <state>\nExample: " << argv[0] << " Toronto 0" << endl;
        return 1;
    }
    int edge_type = 0;
    int state = atoi(argv[2]);
    string suffix;
    if(state == 0) suffix = "MonoGSTBasic";
    else if(state == 6) suffix = "wosusp";
    else if(state == 5) suffix = "womono";
    else if(state == 3) suffix = "tvedge";
    else if(state == 7) suffix = "MonoGSTPlus";
    else suffix = "Other";
    Log::setLogFile("results/" + string(argv[1]) + "_ablation_" + suffix + "_result.txt");
    Log::setConsoleLevel(LogLevel::LOG_INFO);   
    Log::setFileLevel(LogLevel::LOG_IMPORTANT);   

    string graphname = argv[1];
    string dirname = "./data/" + graphname + "/";
    fs_filesystem = dirname;


    auto Go_I2s = [&]<typename edgetype>()
    {
        Graph<edgetype> g;
        auto query = read_query_file();
        Log::debug("Loaded everything");
        Improved2star_Abl<edgetype> i2s;
        i2s.set_state(state);

        for(int i = 0; i < (int)query.size(); i++)
        {
            //if(i!=9) continue;
            int gsize = query[i].size();
            map <int, vector <vector<int>>> query_new;
            for(int j = 0; j < gsize; j++)
            {
                for(auto x: query[i][j])
                {
                    int rt = g.Find(x);
                    if(query_new.find(rt) == query_new.end())
                        query_new[rt] = vector <vector<int>>(gsize, vector<int>());
                    query_new[rt][j].push_back(x);
                }
            }

            Log::info("Query " + to_string(i) + ": ");
            Timer::start("Query " + to_string(i));

            Tree<edgetype> t(numeric_limits<edgetype>::max());
            for(auto [_, qn]: query_new)
            {
                auto tn = i2s.solve(g, qn);
                if(tn.get_sum_weight() < t.get_sum_weight())
                    t = tn;
            }


            auto duration = Timer::stop("Query " + to_string(i), LogLevel::LOG_INFO, false);
            edgetype sum_weight = t.get_sum_weight();
            if(sum_weight == numeric_limits<double>::max()) sum_weight = -1;
            Log::log(LogLevel::LOG_INFO, "Time: " + to_string(duration) + "s" + " sum_weight: " + to_string(sum_weight));
            Log::log(LogLevel::LOG_IMPORTANT, to_string(duration) + " " + to_string(sum_weight));
        }
    };

    if(edge_type == 0) Go_I2s.template operator()<double>();
    else Go_I2s.template operator()<int>();
    
    
    return 0;
} 