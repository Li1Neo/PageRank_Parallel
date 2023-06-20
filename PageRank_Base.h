#ifndef PAGERANK_BASE_H_INCLUDED
#define PAGERANK_BASE_H_INCLUDED
// PageRank 顺序实现

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
using namespace std;


// edges[i] 表示节点 i 指向的所有节点的列表[vector]（出边列表）
// incoming_edges[j] 表示指向节点 j 的所有节点的[vector]（入边列表）
// out_degrees[i] 表示节点 i 的出度[int]
// scores[i] 表示节点i的PageRank分数
vector<double> pagerank_base(int n, const vector<vector<int>>& edges, double damping = 0.85, double tol = 1e-5, int max_iter = 100) {
    // double dampingFactor = 0.85 阻尼系数
    // double tolerance = 1e-5
    vector<double> scores(n, 1.0 / n); // 先把每个节点的PageRank分数初始化为1/n
    vector<int> out_degrees(n);
    vector<vector<int>> incoming_edges(n);
    for (int i = 0; i < n; ++i) {
        out_degrees[i] = edges[i].size(); // 初始化每个顶点的出度
        for (int j : edges[i]) {
            incoming_edges[j].push_back(i); // 初始化incoming_edges
        }
    }
    for (int iter = 0; iter < max_iter; ++iter) { // 迭代max_iter次
        vector<double> prev_scores = scores;
        for (int i = 0; i < n; ++i) {
            double incoming_score = std::accumulate(
                incoming_edges[i].begin(), incoming_edges[i].end(), 0.0,
                [&](double s, int j) { return s + prev_scores[j] / out_degrees[j]; }
            );
            scores[i] = (1 - damping) / n + damping * incoming_score; // 用于解决死胡同问题
        }
        double err = std::inner_product(scores.begin(), scores.end(), prev_scores.begin(), 0.0,
                                [](double a, double b) { return a + b; },
                                [](double a, double b) { return std::abs(a - b); });
        if (err < tol) {
            cout << "iter: "<< iter << ", err: "<< err << endl;
            break;
        }
    }
    return scores;
}

vector<double> pagerank_base2(int n, const vector<vector<int>>& edges, double damping = 0.85, double tol = 1e-5, int max_iter = 1000) {
    // double dampingFactor = 0.85 阻尼系数
    // double tolerance = 1e-5
    vector<double> scores(n, 1.0 / n); // 先把每个节点的PageRank分数初始化为1/n
    vector<double> newpr(n, (1.0 - damping) / n);
    vector<int> out_degrees(n);
    double err = 1.0;
    for(int i = 0; i < n; ++i) {
        out_degrees[i] = edges[i].size(); // 初始化每个顶点的出度
    }
    for(int iter = 0; iter < max_iter; ++iter) { // 迭代max_iter次
        for(int i = 0; i < n; ++i) {
            for(int it=0; it<edges[i].size(); ++it) {
                newpr[edges[i][it]] += damping * (scores[i] / out_degrees[i]);
            }
        }
        err = 0;
        for(int i=0; i<n; ++i) {
            err += fabs(scores[i] - newpr[i]);
			scores[i] = newpr[i];
			newpr[i] = (1.0 - damping) / n;
        }
        if(err < tol) {
            cout << "iter: "<< iter << ", err: "<< err << endl;
            break;
        }
    }
    return scores;
}

#endif