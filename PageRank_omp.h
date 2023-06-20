/*** 
 * @Description  : Pagerank OpenMP实现
 * @Version      : 1.0
 * @Author       : LilNeo
 * @Date         : 2023-05-26 15:08:25
 * @LastEditors  : wy
 * @LastEditTime : 2023-05-28 22:24:05
 * @FilePath     : /Pagerank/PageRank_omp.h
 * @Copyright 2023 LilNeo, All Rights Reserved. 
 */
#ifndef PAGERANK_OMP_H_INCLUDED
#define PAGERANK_OMP_H_INCLUDED

#include <omp.h>
// #include <cmath>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#define thread_num 8
#define chunk 2
using std::vector, std::cout, std::endl;

vector<double> pagerank_omp(int n, const vector<vector<int>>& edges, double damping = 0.85, double tol = 1e-5, int max_iter = 1000) {
    vector<double> scores(n, 1.0 / n);
    vector<int> out_degrees(n);
    vector<double> newpr(n, (1.0 - damping) / n);
    double err = 1.0;
    #pragma omp parallel num_threads(thread_num)
    {
        // #pragma omp for schedule(dynamic, chunk)
        #pragma omp for schedule(guided)
        for(int i = 0; i < n; ++i) {
            out_degrees[i] = edges[i].size();
        }
        while(err > tol) {
            // #pragma omp for schedule(dynamic, chunk)
            #pragma omp for schedule(guided)
            for (int i = 0; i < n; ++i) {
                for (int it = 0; it < edges[i].size(); ++it) {
                    #pragma omp atomic
                    newpr[edges[i][it]] += damping * (scores[i] / out_degrees[i]);
                }
            }
            err = 0;
            #pragma omp for reduction(+:err)
            for (int i = 0; i < n; ++i) {
                err += fabs(scores[i] - newpr[i]);
                scores[i] = newpr[i];
                newpr[i] = (1.0 - damping) / n;
            }
        }
    }
    cout << "err:" << err << endl;
    return scores;
}

vector<double> pagerank_omp2(int n, const vector<vector<int>>& edges, double damping = 0.85, double tol = 1e-5, int max_iter = 100) {
    // double dampingFactor = 0.85 阻尼系数
    // double tolerance = 1e-5
    vector<double> scores(n, 1.0 / n); // 先把每个节点的PageRank分数初始化为1/n
    vector<int> out_degrees(n);
    vector<vector<int>> incoming_edges(n);
    double err = 1.0;
    // #pragma omp parallel num_threads(thread_num) shared(incoming_edges)
    {
        // #pragma omp single
        for (int i = 0; i < n; ++i) {
            out_degrees[i] = edges[i].size(); // 初始化每个顶点的出度
            for(int j : edges[i]) {
                incoming_edges[j].push_back(i); // 初始化incoming_edges
            }
        }
        #pragma omp parallel num_threads(thread_num)
        while(err > tol) {
            vector<double> prev_scores = scores;
            #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                double incoming_score = 0.0;
                for (int j : incoming_edges[i]) {
                    incoming_score += prev_scores[j] / out_degrees[j];
                }
                scores[i] = (1 - damping) / n + damping * incoming_score; // 用于解决死胡同问题
            }
            err = std::inner_product(scores.begin(), scores.end(), prev_scores.begin(), 0.0,
                                    [](double a, double b) { return a + b; },
                                    [](double a, double b) { return std::abs(a - b); });
        }
    }
    return scores;
}


std::vector<double> pagerank(int n, const std::vector<std::vector<int>>& edges, double damping = 0.85,
                            double tol = 1e-5, int max_iter = 100) {
    std::vector<double> scores(n, 1.0 / n);
    std::vector<int> out_degrees(n);
    std::vector<std::vector<int>> incoming_edges(n);
    double err = 1.0;

    // #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        out_degrees[i] = edges[i].size();
        for (int j : edges[i]) {
            incoming_edges[j].push_back(i);
        }
    }

    while (err > tol && max_iter > 0) {
        std::vector<double> prev_scores = scores;

        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double incoming_score = 0.0;
            for (int j : incoming_edges[i]) {
                incoming_score += prev_scores[j] / out_degrees[j];
            }
            scores[i] = (1 - damping) / n + damping * incoming_score;
        }

        err = std::inner_product(scores.begin(), scores.end(), prev_scores.begin(), 0.0,
                                [](double a, double b) { return a + b; },
                                [](double a, double b) { return std::abs(a - b); });

        --max_iter;
    }

    return scores;
}


#endif