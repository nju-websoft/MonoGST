#ifndef LOADQUERY_H
#define LOADQUERY_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <map>
#include <random>

#include "AdjacencyList.hpp"
//Load for KeyKG
template<typename NodeType>
class LoadQuery
{
    string filepath;
    std::mt19937 gen;
    
    public:
        vector <vector <vector <NodeType> > > Queryset;
        int seed;
        NodeType GenerateRandomNumber(NodeType l, NodeType r)
        {
            std::uniform_int_distribution<NodeType> rnd(l, r);
            return rnd(gen);
        }
        // g: 1~10 small: gl,gr: 2 10  mid: gl,gr: 20 50  large: gl,gr: 100 1000
        void GenerateRandomQuery(int TestCase, NodeType g, NodeType gl, NodeType gr, NodeType n, int seed = rand())
        {
            this->seed = seed;
            gen.seed(seed);
            for(int i = 0; i < TestCase; i++)
            {
                vector <vector <NodeType> > Query;
                for(NodeType j = 0; j < g; j++)
                {
                    vector <NodeType> v;
                    NodeType cnt = GenerateRandomNumber(gl, gr);
                    for(NodeType k = 0; k < cnt; k++)
                        v.push_back(GenerateRandomNumber(0, n-1));
                    Query.push_back(v);
                }
                Queryset.push_back(Query);
            }
        }
        
        LoadQuery(string fp, int type): filepath(fp)
        {
            if(type == 0)
            {
                string filename = filepath + "/queryList.txt";
                ifstream inputFile(filename);
                if (!inputFile.is_open()) 
                    cerr << "Can't open file: " << filename << endl;
                string line;
                while (getline(inputFile, line)) 
                {
                    vector <vector <int> > Query;
                    istringstream lines(line);
                    string lin;
                    while(getline(lines, lin, ','))
                    {
                        auto num = extractIntegers(lin);
                        Query.push_back(num);
                    }    
                    Queryset.push_back(Query);
                }
                inputFile.close();
            }
            else if(type == 1)
            {
                auto knm_kid = get_knm_to_kid();
                auto kid_vid = get_kid_to_vid();
                string filename = filepath + "/query.txt";
                ifstream inputFile(filename);
                if (!inputFile.is_open()) 
                    cerr << "Can't open query: " << filename << endl;
                string line;
                int testcase = 0;
                while (getline(inputFile, line)) 
                {
                    vector <vector <int> > Query;
                    istringstream iss(line);
                    string knm;
                    while(iss >> knm)
                    {
                        knm = transl(knm);
                        if(knm_kid.find(knm) == knm_kid.end())
                        {
                            cerr << "Name " << knm << " don't hava kid!" << endl;
                            continue;
                        }
                        if(kid_vid.find(knm_kid[knm]) == kid_vid.end())
                        {
                            cerr << "Kid " << knm_kid[knm] << "don't have node set!" << endl;
                            continue;
                        }
                        Query.push_back(kid_vid[knm_kid[knm]]);
                    }
                    int sumsize = 0;
                    for(int i = 0; i < Query.size(); i++)
                        sumsize += Query[i].size();
                    cerr << "Average size of query " << testcase << ": " << sumsize / (1.0*Query.size()) << "\n";
                    testcase++;
                    Queryset.push_back(Query);
                }
                inputFile.close();
            }
        }
    
    private:
        std::vector<NodeType> extractIntegers(const std::string& input) 
        {
            std::vector<NodeType> integers;
            std::istringstream iss(input);
            std::string token;

            while (iss >> token) 
            {
                try 
                {
                    NodeType number = stoi(token);
                    integers.push_back(number);
                } catch (const std::invalid_argument& e) 
                {
                    std::cerr << "not a number" << std::endl;
                }
            }

            return integers;
        }


        /*
        real
        nodeName.txt: Mapping from vertex ID to vertex name (i.e., entity URI).
        query.txt: Each line is a keyword query containing a set of keyword names.
        kwName.txt: Mapping from keyword ID to keyword name.
        kwMap.txt: Mapping from keyword ID to vertex IDs. The first value of each line is keyword ID, and the rest are vertex IDs.
        keyword ID -> keyword name
        keyword ID -> vertex IDs
        */
       
        string transl(string str)
        {
            transform(str.begin(), str.end(), str.begin(), [](unsigned char c){return tolower(c);});
            return str;
        }
        map <string, int> get_knm_to_kid()
        {
            string filename = filepath + "/kwName.txt";
            ifstream inputFile(filename);
            if (!inputFile.is_open()) 
                cerr << "Can't open kwName: " << filename << endl;
            
            map <string, int> knm_kid;
            string line;
            while (getline(inputFile, line)) 
            {
                //cerr << line << endl;
                istringstream iss(line);
                int kid;
                string knm;
                iss >> kid >> knm;   
                knm = transl(knm);
                knm_kid[knm] = kid;
            }
            inputFile.close();
            return knm_kid;
        }

        map <int, vector <int> > get_kid_to_vid()
        {
            string filename = filepath + "/kwMap.txt";
            ifstream inputFile(filename);
            if (!inputFile.is_open()) 
                cerr << "Can't open kwMap: " << filename << endl;
            
            map <int, vector <int> > kid_vid;
            string line;
            while (getline(inputFile, line)) 
            {
                istringstream iss(line);
                int kid, vid;
                iss >> kid;
                vector <int> svid;
                while(iss >> vid)
                    svid.push_back(vid);
                kid_vid[kid] = svid;
            }
            inputFile.close();
            return kid_vid;
        }
};

#endif