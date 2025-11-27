#ifndef UTILITIES_H
#define UTILITIES_H

#include <cassert>
#include <iostream>
#include <deque>
#include <limits>
#include <random>
#include <algorithm>
#include <unordered_set>

#include "AdjacencyList.hpp"

/*
N: vextex number. 
M: edge number. 
w: edge range. if DistanceType is int and w = 1 -> unweighed graph
perc: random choose vextex f(N) = f(N*perc) for p=1/2 or random choose v in [N*perc, N-1] with equal probability for p=1/2
seed: random seed
*/
template<typename DistanceType>
inline AdjacencyList<int, DistanceType> GenGraph(int N, int M, DistanceType w, double perc, int seed)
{
    std::mt19937 gen(seed);

    std::uniform_int_distribution<int> rnd(0, std::numeric_limits<int>::max());
    std::uniform_real_distribution<double> rndw(0, 1);

    AdjacencyList<int, DistanceType> adjList(N);

    auto randvex = [&](auto &&randvex, int x) -> int
    {
        int y=x*perc;
        if(x <= 6 || y == 0 || y == x) return rnd(gen)%x;
        if(rndw(gen)<=0.5) return rnd(gen)%(x-y) + y;//[y,x-1]
        else return randvex(randvex, y);// [0, y-1]
    };
    auto gene_edge_value = [&]()
    {
        DistanceType r;
        if(std::is_same<DistanceType, int>::value) 
        {
            int rw=w;
            if(w==1) r = 1;
            else r = rnd(gen)%(rw+1)+1;
        } 
        else r = rndw(gen)*w;
        return r;
    };

    std::vector <int> permu;
    for(int i = 0; i < N; i++) permu.push_back(i);
    shuffle(permu.begin(), permu.end(), gen);
    for(int i = 0; i < M; i++)
    {
        int u = randvex(randvex, N), v = randvex(randvex, N);
        u = permu[u], v = permu[v];
        adjList.AddUndirectedEdge(u, v, gene_edge_value());
    }

    return adjList;
}

#endif