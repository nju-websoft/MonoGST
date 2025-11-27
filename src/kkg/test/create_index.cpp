#include "../HL/HL_all.hpp"
#include "../AdjacencyList.hpp"
#include "../Utilities.hpp"
#include "../DistanceOracle.hpp"
#include "../Load_Graph.hpp"
#include "../Timer.hpp"
#include "../Prepare_Test.hpp"

#include <type_traits>
#include <format>
#include <iostream>
#include <sstream>
#include <iomanip>

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

    string graph_file = "./data/" + graphname;

    FileManager fop("KeyKG_Index/" + graphname);
    fop.clearFile("result.txt");

    auto GoIndex = [&]<typename EdgeType>()
    {
        LoadGraph<int, EdgeType> g(graph_file, graph_type == "0" ? 0 : 1, 1);

        DistanceOracle<int, EdgeType, HL_all> DO(g.adjList, 1);
        fop.writeToFile("result.txt", "IndexTime: " + to_string(DO.oracle.cost_time) + " ms.");
        DO.oracle.Output_Index("KeyKG_Index/" + graphname + "/");
    };
    
    if(graph_type == "0") 
        GoIndex.template operator()<int>();
    else
        GoIndex.template operator()<long double>();

    return 0;
}