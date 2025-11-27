#ifndef Prepare_Test_H
#define Prepare_Test_H

#include <cassert>
#include <iostream>
#include <limits>
#include <random>
#include <chrono>
#include <algorithm>
#include <variant>
#include <fstream>
#include <filesystem>
#include <map>
#include <set>
using std::set;
using std::pair; 
using std::vector;
using std::map;
using std::string;
using std::cerr;
using std::to_string;
using std::tuple;
using Clock = std::chrono::high_resolution_clock;
namespace fs = std::filesystem;
class FileManager 
{
private:
    std::string filepath;

    void createDirectoryIfNotExist() 
    {
        if (!fs::exists(filepath)) 
        {
            std::cout << "Directory does not exist. Creating it..." << std::endl;
            fs::create_directories(filepath);
        }
    }

public:
    FileManager(const std::string& path) : filepath(path) 
    {
        createDirectoryIfNotExist();
    }

    void writeToFile(const std::string& filename, const std::string& info) 
    {
        std::string file = filepath + "/" + filename;

        std::ofstream outFile(file, std::ios::app);
        if (!outFile) 
        {
            std::cerr << "Failed to open the file." << std::endl;
            return;
        }

        outFile << info << '\n';
    }

    void writeToFile(const std::string& filename, const vector <string>& Info) 
    {
        std::string file = filepath + "/" + filename;

        std::ofstream outFile(file, std::ios::app);
        if (!outFile) 
        {
            std::cerr << "Failed to open the file." << std::endl;
            return;
        }
        for(auto info: Info)
            outFile << info << '\n';
    }

    void clearFile(const std::string& filename) {
        std::string file = filepath + "/" + filename;
        std::ofstream outFile(file, std::ios::trunc); 
    }
};
//NY BAY COL FLA NW NE NW CAL LKS E W CTR USA
class Parameter
{
public:
    map <string, int> graph_offset;
    map <string, string> graph_type;
    set <string> DIMACS_S;
    map <string, int> g_query_type;
    
    void Load_Parameters()
    {
        vector <string> DIMACS_List;
        DIMACS_List.push_back("NY");
        DIMACS_List.push_back("BAY");
        DIMACS_List.push_back("COL");
        DIMACS_List.push_back("FLA");
        DIMACS_List.push_back("NW");
        DIMACS_List.push_back("NE");
        DIMACS_List.push_back("CAL");
        DIMACS_List.push_back("LKS");
        DIMACS_List.push_back("E");
        DIMACS_List.push_back("W");
        DIMACS_List.push_back("CTR");
        DIMACS_List.push_back("USA");
        for(auto v: DIMACS_List)
        {
            graph_offset[v] = 1;
            DIMACS_S.insert(v);
        }

        vector <string> graph_list;
        graph_list.push_back("Toronto");
        graph_list.push_back("Movielens");
        graph_list.push_back("DBLP");
        graph_list.push_back("Dbpedia");
        graph_list.push_back("LinkedMDB");
        graph_list.push_back("LUBM-50K");
        graph_list.push_back("LUBM-500K");
        graph_list.push_back("LUBM-5M");
        graph_list.push_back("Mondial");
        graph_list.push_back("OpenCyc");
        graph_list.push_back("YAGO");
        graph_list.push_back("example");

        graph_list.push_back("baidu");
        graph_list.push_back("epinions");

        for(auto gn: graph_list) graph_offset[gn] = 0;
        graph_offset["Toronto"] = 1;
        graph_offset["Movielens"] = 1;
        graph_offset["DBLP"] = 1;

        for(auto gn: graph_list) g_query_type[gn] = 1;
        g_query_type["LUBM-50K"] = 0;
        g_query_type["LUBM-500K"] = 0;
        g_query_type["LUBM-5M"] = 0;
        g_query_type["example"] = 0;

        graph_type["0"] = graph_type["UW"] = "UW";
        graph_type["1"] = graph_type["IW"] = "IW";
        graph_type["TG"] = "Travel_Time";
        graph_type["DG"] = "Distance";
    }
};


#endif