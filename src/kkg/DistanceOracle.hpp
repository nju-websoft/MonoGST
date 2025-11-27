#ifndef DISTANCE_ORACLE_HPP
#define DISTANCE_ORACLE_HPP
#include "AdjacencyList.hpp"
#include "Utilities.hpp"
using std::vector;

template<typename NodeType, typename DistanceType, template<typename, typename> class OracleType>
class DistanceOracle
{
public:
    using AdjList = AdjacencyList<NodeType, DistanceType>;
    OracleType<NodeType, DistanceType> oracle; 

    DistanceOracle(AdjList& lis): oracle(lis) {}

    DistanceOracle(AdjList& lis, string filepre): oracle(lis, filepre) {}

    DistanceOracle(AdjList& lis, string filepre, int is_q_p): oracle(lis, filepre, is_q_p) {}

    DistanceOracle(AdjList& lis, string filepre, int is_q_p, int is_read_index, int is_ternary_search, int need_random_shuffle): oracle(lis, filepre, is_q_p, is_read_index, is_ternary_search, need_random_shuffle) {}
    
    DistanceOracle(AdjList& lis, int K): oracle(lis, K){}

    DistanceOracle(AdjList& lis, int K, string filepre): oracle(lis, K, filepre) {}

    DistanceOracle(AdjList& lis, int K, int is_q_p): oracle(lis, K, is_q_p){}

    DistanceOracle(AdjList& lis, double part_one, int is_q_p): oracle(lis, part_one, is_q_p) {}

    DistanceOracle(AdjList& lis, string filepre, double part_one, int is_q_p): oracle(lis, filepre, part_one, is_q_p) {}

    DistanceOracle(AdjList& lis, int K, string filepre, int is_q_p): oracle(lis, K, filepre, is_q_p) {}

    DistanceType Query_Distance(NodeType u, NodeType v) 
    {
        return oracle.Query(u, v);
    }

    vector<NodeType> Query_Path(NodeType u, NodeType v) 
    {
        return oracle.Query_Path(u, v);
    }
};

#endif  // DISTANCE_ORACLE_HPP