# A Practical Sublinear Approximation for Group Steiner Tree

This is the source code of the paper 'A Practical Sublinear Approximation for Group Steiner Tree'.

## Environment

C++ (need to support C++20)

## Data

Our data is available from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18425871.svg)](https://doi.org/10.5281/zenodo.18425871).

Extract all the .zip files into the directory `data`, you could see `data/example` for an example. 

Each data contains a `graph.txt` and a `query.txt` used in our experiments, including 5 real-world datasets `Toronto`, `MovieLens`, `DBLP`, `LinkedMDB` and `DBpedia`.

Each data directory contains 2 files, including:

- `graph.txt`: the first line contains two values 'n','m', which is the number of vertices and the number of edges in the graph. Then next m lines contains three values 'u', 'v', 'w' which means there is an undirected edge between 'u' and 'v' weighted by 'w' which is computed by the Informativeness-based Weighting (IW) scheme.
- `query.txt`: The first value is the number of queries. For each query, the first value is the number of groups 'g'. For the next g lines, the first value is the size 'f' of the group, then the next f values are the vertices in this group. For example, the file `data/example/query.txt` means there is 1 query; in this query, g = 4 groups. The 4 groups are: {1}, {5, 7}, {10, 11}, and {4, 6, 8}. 

## Compile Command

- `make MonoGST`  or `g++ -std=c++20 -O2 -g -Isrc/MonoGST/include -DDEBUG -o MonoGST src/MonoGST/tests/test_Improved2star.cpp src/MonoGST/src/GlobalUtils.cpp src/MonoGST/src/Log.cpp` to compile MonoGST and MonoGST+.
- `make ImprovAPP` or `g++ -std=c++20 -O2 -g -o ImprovAPP src/ImprovAPP/ImprovAPP.cpp` to compile the ImprovAPP.
- `make PrunedDP` or `g++ -std=c++20 -O2 -g -o PrunedDP src/PrunedDP/PrunedDP.cpp` to compile the PrunedDP.
- `make PartialOPT` or `g++ -std=c++20 -O2 -g -o PartialOPT src/PartialOPT/PartialOPT.cpp` to compile the PartialOPT.
- `make KKG_Index` or `g++ -std=c++20 -O2 -g -o KKG_Index src/kkg/test/create_index.cpp` to compile the Index-generator for KeyKG+; `make KKG_Run` or `g++ -std=c++20 -O2 -g -o KKG_Run src/kkg/test/run_kkg.cpp` to compile the KeyKG+ itself.
- `make 2starh` or `g++ -std=c++20 -O2 -g -o 2starh src/2starh/2starh.cpp` to compile the 2-star heuristic.
- `make DST` to compile the DST. Currently, we have only conducted experiments with this algorithm in an Ubuntu environment. According to the requirements specified in the original author's paper, the necessary environment is as follows: C++ 17 or later, Boost, CMake. In Ubuntu, these can be installed via
    ```
    apt install -y make \
                build-essential \
                cmake \
                libboost-all-dev
    ```
    

### Run Algorithms

The data format and result format are identical in all scenarios described in the paper, and we take `data/example/` as an example.

Run MonoGST+ for `example`: `./MonoGST example`.

Run MonoGST for `example`: `./MonoGST example 0` .

Run the ablation experiments of MonoGST+ for `example`: `./MonoGST example x`, `x` is an integer parameter to control the variants of MonoGST+. If x=6, it will run `w/o suspension`; if x=4, it will run `w/o suspension & filtering`; if x=3, it will run `w/o re-selection`; if x=7, it will run `MonoGST+`.

Before run KeyKG+, you need first build the index of the graph, use `./KKG_Index example` to finish it, then you could use `./KKG_Run example` to answer all queries for the graph.

Run DST for `example`: `./src/DST/build/main/Main example`.

Run other baseline `b` for `example`: `./b example`, you need replace "`b`" to the name of the baseline of  ` ImprovAPP`, `PrunedDP`, `PartialOPT` or `2starh`.

Run other baselines `b` for `example`: `./b example`. Replace "`b`" with the name of the baseline: `ImprovAPP`, `PrunedDP`, `PartialOPT`, or `2starh`.

The result of the algorithm `b` for example will be replaced in directory `results/example_b_result.txt`, each line of the text outputs the result of an inquiry, with the first value being the time (in seconds) and the second value being the sum of the edge weights in GST. Specifically, the index of the KeyKG+ and its related information will be output to the directory `KeyKG_index/example/index.bin`.

## Citation

If you think our algorithms or our experimental results are useful, please kindly cite our paper.