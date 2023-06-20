// PageRank,基于文件块更新的方式
// 指定数据以长度为N进行分块，存取
#include <cassert>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<unordered_map>
#include<unordered_set>
#include<iomanip>
#include<time.h>
#include<math.h>
#include<stdio.h>
#include "helper_utils.h"
using namespace std;


// 写入块
void WriteBlock(int num_vex, const int idx, vector<int>& linked_vex, int block_size) {
	vector<vector<int>> blocks(num_vex/block_size + 1);
	unordered_set<int> to_write_blockno_list;
	for(auto i : linked_vex) {
		int blockno = i / block_size;
		blocks[blockno].push_back(i); // 第blockno个block存放的节点列表
		if(to_write_blockno_list.find(blockno) == to_write_blockno_list.end()) {
			to_write_blockno_list.insert(blockno);
		}
	}
	ofstream out;
	char block_file[20];
	for(int i : to_write_blockno_list) {
		sprintf(block_file, "block/block-%d.txt", i);
		out.open(block_file, ios::app);
		out << idx << "\t" << linked_vex.size() << "\t" << blocks[i].size();
		for(auto j : blocks[i]) {
			out << "\t" << j;
		}
		out.put('\n');
		out.close();
	}
}

// 读取原始数据构造Block
void InitBlock(int num_vex, const string path, int block_size) {
	ifstream infile;
	infile.open(path);
	int u, v;
	// 建图
	vector<vector<int>> edges(num_vex);
	while(infile>>u>>v) {
		edges[u].push_back(v);
	}
	infile.close();
	// 对每个edges[i]sort一下，让位于一个stripe的节点相邻，减少io次数
	for(int i=0; i<num_vex; ++i) {
		sort(edges[i].begin(), edges[i].end());
		WriteBlock(num_vex, i, edges[i], block_size);
	}
	edges.clear();
}

vector<double> PageRank_Block_Stripe(int n, const string path, unordered_map<int, int> &idx2id, int block_size = 100, double damping = 0.85, double tol = 1e-6, int max_iter = 1000) {
    InitBlock(n, path, block_size);
    vector<double> newpr(block_size, (1.0 - damping) / n);
    vector<double> scores(n, 1.0 / n);
    vector<int> out_degrees(n);
    double err = 1.0;
    int degree, size, cnt = 0;
    ifstream block_infile;
    ofstream newpr_outfile;
    ifstream newpr_infile;
    char block_file[20];
    int u, v, num_blocks = n % block_size == 0 ? (n / block_size) : (n / block_size + 1);
    for(int iter=0; iter<max_iter; ++iter) {
        err = 0;
        newpr_outfile.open("newpr.txt", ios::trunc);
        newpr_outfile.close();
        newpr_outfile.open("newpr.txt", ios::app);
        for(int i = 0; i < num_blocks; ++i) { //对每个块分别更新
            // i 为块号
            sprintf(block_file, "block/block-%d.txt", i);
            block_infile.open(block_file);
            while(block_infile >> u >> degree >> size) {
                while(size--) {
                    block_infile >> v;
                    assert(v/block_size == i);
                    newpr[v % block_size] += damping * scores[u] / degree;
                }
            }
            block_infile.close();
            for(int j=0; j<block_size; ++j) {
                newpr_outfile << newpr[j] << endl;
                newpr[j] = (1.0 - damping) / n;
            }
        }
        ++cnt;
        newpr_infile.open("newpr.txt");
        for(int i = 0; i < n; i++) {
            newpr_infile >> newpr[i%block_size];
            int blockno = i/block_size;
            err += fabs(newpr[i%block_size] - scores[i]);
            scores[i] = newpr[i%block_size];
        }
        for(int i = 0; i < block_size; ++i) {
            newpr[i] = (1.0 - damping) / n;
        }
        if(err < tol) {
            cout<<"迭代次数:" << cnt <<",err:"<<err<<endl;
            break;
        }
        newpr_outfile.close();
        newpr_infile.close();
    }
    return scores;
}

int main(int argc, char *argv[]) {
	// system("mkdir block");
	system("rm block/*.txt");
	const string path = "newData.txt";
	int n = 6263; // 节点总数
    int block_size = 10;
    double damping = 0.85;
    double tol = 1e-6;
    cout<<"分块大小："<<block_size<<"，damping:" <<damping <<"，tol:"<<tol<<endl;
	unordered_map<int, int> idx2id;
	ifstream infile("idx2id.txt");
    if (!infile) {
        cout << "Error opening file." << endl;
        exit(-1);
    }
    int key;
    int value;
    while(infile >> key >> value) {
        idx2id[key] = value;
    }
    infile.close();
	clock_t start, end;

    start = clock();
    vector<double> scores = PageRank_Block_Stripe(n, path, idx2id, block_size=block_size, damping=damping, tol=tol);
    end = clock();
    float total_time = (end - start) / float(CLOCKS_PER_SEC);
	cout << "Time: " << total_time << endl;
	vector<pair<int, double> > indexed_scores;
    for(int i=0; i<scores.size(); ++i) {
        indexed_scores.push_back({i, scores[i]});
    }
    sort(indexed_scores.begin(), indexed_scores.end(), cmp);
    print_pagerank_scores(indexed_scores, idx2id, 100);
    return 0;
}