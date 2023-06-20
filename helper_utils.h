#ifndef HELPER_UTILS_H_INCLUDED
#define HELPER_UTILS_H_INCLUDED

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
using std::vector, std::string, std::ifstream, std::ofstream, std::cout, std::endl, std::unordered_map, std::pair;

// 根据Data.txt建图
//vector<vector<int>> build_edges(string path) {
//    int u, v;
//    unordered_map<int, int> um;
//    vector<vector<int>> edges;
//    ifstream infile;
//    infile.open(path);
//    if (!infile) {
//        cout << "无法打开"<<path<<"!" << endl;
//        exit(1);
//    }
//    int i = 0;
//    while(infile>>u>>v) {
//        cout<<u<<" "<<v<<endl;
//        if(um.find(u) != map.end()) {
//            um.insert({u, i});
//            ++i;
//            edges.push_back(vector<int>());
//        }
//        if(um.find(v) != map.end()) {
//            um.insert({v, i});
//            ++i;
//        }
//        edges.push_back(vector<int>());
//        edges[um[u]].push_back(v);
//        hashtable.push_back(u);
//
//        } else {
//            edges[u].push_back(v);
//        }
//        cout<<endl;
//        cout<<edges[0].size()<<endl;
//        cout<<edges[1].size()<<endl;
//        cout<<edges[u].size()<<endl;
//        cout<<edges.size()<<endl;
//    }
//    infile.close();//关闭文件
//    return edges;
//}

vector<vector<int>> build_edges(int num_vex, unordered_map<int, int> &um, string path) {
    int u, v;
    vector<vector<int>> edges(num_vex);
    ifstream infile;
    ofstream outfile;
    outfile.open("newData.txt");
    infile.open(path);
    if (!infile) {
        cout << "无法打开"<<path<<"!" << endl;
        exit(1);
    }
    int i = 0;
    while(infile>>u>>v) {
        outfile << um[u] << ' ' << um[v] << endl;
        edges[um[u]].push_back(um[v]);
        // edges[um[v]].push_back(um[u]);
    }
    infile.close();//关闭文件
    outfile.close();
    return edges;
}

// Print the PageRank scores for each node
void print_pagerank_scores(const vector<pair<int, double>>& indexed_scores, unordered_map<int, int> &idx2id, int topk) {
    cout << "PageRank scores(top" << topk << "):" << endl;
    for (int i = 0; i < topk; i++) {
        cout << "Node " << idx2id[indexed_scores[i].first] << ": " << indexed_scores[i].second << std::endl;
    }
}


bool cmp(const pair<int, double>& a, const pair<int, double>& b) {
    return a.second > b.second;
}




#endif // HELPER_UTILS_H_INCLUDED
