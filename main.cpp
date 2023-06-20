

/*** 
 * @Description  : To be filled
 * @Version      : 1.0
 * @Author       : LilNeo
 * @Date         : 2023-04-21 23:03:49
 * @LastEditors  : wy
 * @LastEditTime : 2023-05-28 17:08:33
 * @FilePath     : /Pagerank/main.cpp
 * @Copyright 2023 LilNeo, All Rights Reserved. 
 */
#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include "PageRank_Base.h"
#include "PageRank_omp.h"
#include "helper_utils.h"
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <omp.h>
using std::vector;
using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
    ifstream infile;
    int u, v;
    unordered_map<int, int> um;
    unordered_map<int, int> idx2id;
    string path = argv[1]; // "web-Google.txt"或"Data.txt"
    cout << "Processing " << path << "..." << endl;
    infile.open(path);
    int maxn = 0;
    int minn = 0x7fffffff;
    int num_edge = 0;
    int num_vex = 0;
    int cnt = 0;
    while(infile>>u>>v) {
        if(um.find(u)==um.end()) {
            um[u] = num_vex++;
        }
        if(um.find(v)==um.end()) {
            um[v] = num_vex++;
        }
        ++num_edge;
        maxn = max(max(maxn, u),v);
        minn = min(min(minn, u),v);
    }    
    for(auto kv : um) {
        idx2id[kv.second] = kv.first;
    }
    ofstream outfile("idx2id.txt");
    if (!outfile) {
        cout << "Error opening file." << endl;
        exit(-1);
    }
    for (auto it = idx2id.begin(); it != idx2id.end(); it++) {
        outfile << it->first << " " << it->second << endl;
    }
    outfile.close();

    cout<<"最大节点为："<<maxn<<endl;
    cout<<"最小节点为："<<minn<<endl;
    cout<<"节点个数："<<num_vex<<"，边个数："<<num_edge<<endl;
    infile.close();
    vector<vector<int>> edges = build_edges(num_vex ,um, path);
    // cout<<edges.size()<<endl;
    // for(int i=0; i<edges.size(); ++i) {
    //     cout<<edges[i].size()<<'\t';
    //     for(int j=0; j<edges[i].size(); ++j) {
    //         cout<<edges[i][j]<<' ';
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    double damping = 0.85;
    double tol = 1e-6;
    cout<<"damping:" << damping <<"，tol:"<<tol<<endl;
    // clock_t start, finish;
    double seconds;
    // start = clock();
    double start = omp_get_wtime();
    vector<double> scores;
    string opt = argv[2];
    if(opt == "serial") {
        cout << "执行串行程序！" << endl;
        scores = pagerank_base(num_vex, edges, damping=damping, tol=tol);
        // scores = pagerank_base2(num_vex, edges, damping=damping, tol=tol);
    } else if(opt == "omp") {
        cout << "执行OpenMP程序！" << endl;
        scores = pagerank_omp(num_vex, edges, damping=damping, tol=tol);
        // scores = pagerank_omp2(num_vex, edges, damping=damping, tol=tol);
        // scores = pagerank(num_vex, edges, damping=damping, tol=tol);
    } else {
        cout << "未识别的指令！" << endl;
        return -1;
    }
    // finish = clock();
    double finish = omp_get_wtime();
    seconds = finish - start;
    // seconds = (finish - start) / double(CLOCKS_PER_SEC);
    vector<pair<int, double> > indexed_scores;
    for(int i=0; i<scores.size(); ++i) {
        indexed_scores.push_back({i, scores[i]});
    }
    sort(indexed_scores.begin(), indexed_scores.end(), cmp);
    print_pagerank_scores(indexed_scores, idx2id, 10);
    cout<< "用时" << seconds << "s" << endl;
    return 0;
}
/*
// shell
/opt/homebrew/opt/llvm/bin/clang++ -o2 -std=c++17 -fopenmp main.cpp -o main
// Serial + small Dataset:
./main Data.txt serial
// Serial + large Dataset:
./main web-Google.txt serial
// OpenMP + small Dataset:
./main Data.txt omp
// OpenMP + large Dataset:
./main web-Google.txt omp 
*/

/*
最大节点为：8297
最小节点为：3
节点个数：6263，边个数：83852
damping:0.85，tol:1e-06
err:9.66161e-07
PageRank scores:
Node 4037: 0.00219573
Node 2625: 0.00179141
Node 6634: 0.00163986
Node 15: 0.00136354
Node 2398: 0.00123883
Node 2328: 0.00113335
Node 2470: 0.00111939
Node 3089: 0.00108733
Node 6946: 0.00104464
Node 3352: 0.00103998
Node 5412: 0.00103214
Node 4191: 0.00100879
Node 7632: 0.000997653
Node 7553: 0.000959837
Node 737: 0.00094348
Node 1297: 0.000940077
Node 3456: 0.00092863
Node 2237: 0.000926924
Node 5254: 0.00091461
Node 6832: 0.000911348
Node 2066: 0.000880003
Node 4712: 0.00083788
Node 762: 0.000828038
Node 7092: 0.00082663
Node 1186: 0.000820637
Node 4310: 0.000817672
Node 6774: 0.000781911
Node 7620: 0.000776142
Node 2958: 0.000761259
Node 4335: 0.000741297
Node 993: 0.000726703
Node 4828: 0.000724514
Node 3537: 0.000723053
Node 4875: 0.000721428
Node 6006: 0.000718089
Node 2657: 0.000714498
Node 271: 0.000714353
Node 665: 0.00070989
Node 1549: 0.000708851
Node 4735: 0.000700377
Node 4256: 0.000697869
Node 3238: 0.000669153
Node 5484: 0.00065371
Node 825: 0.000653375
Node 3498: 0.000653211
Node 2565: 0.000652006
Node 4261: 0.000646866
Node 3568: 0.000638458
Node 5123: 0.000637808
Node 2654: 0.000631511
Node 3084: 0.000629765
Node 2485: 0.000618246
Node 5079: 0.000615844
Node 6784: 0.000615455
Node 5404: 0.000608415
Node 2535: 0.000608367
Node 2871: 0.000607508
Node 5543: 0.00059846
Node 6334: 0.000596227
Node 5459: 0.000591507
Node 3897: 0.000591286
Node 28: 0.000589809
Node 8042: 0.000579452
Node 4400: 0.000574838
Node 1842: 0.000570951
Node 2859: 0.00056899
Node 6059: 0.000566821
Node 3562: 0.000566263
Node 4600: 0.000564145
Node 2746: 0.000562389
Node 3334: 0.000562291
Node 2576: 0.000562224
Node 3034: 0.000560391
Node 5226: 0.000557964
Node 7809: 0.000555163
Node 2651: 0.000546414
Node 5022: 0.000541817
Node 3321: 0.000525762
Node 1633: 0.000520929
Node 5563: 0.000520435
Node 4099: 0.000517954
Node 1211: 0.000515761
Node 1754: 0.000511208
Node 3459: 0.000509948
Node 8293: 0.000507545
Node 4040: 0.000506732
Node 4666: 0.000505022
Node 1726: 0.000502604
Node 4981: 0.000502497
Node 7961: 0.000501923
Node 4531: 0.000499872
Node 2516: 0.000498897
Node 6124: 0.000497983
Node 6330: 0.000495683
Node 5605: 0.000493932
Node 3962: 0.000492121
Node 3005: 0.000488473
Node 7890: 0.000488144
Node 86: 0.000485411
Node 7214: 0.000483894
用时0.022567s
*/