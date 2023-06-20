#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <mpi.h>
#include "helper_utils.h"
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/mpi.hpp>
using std::vector, std::max, std::min;

namespace boost {
    namespace serialization {
        template<class Archive>
        void serialize(Archive& ar, std::vector<std::vector<int>>& edges, const unsigned int version)
        {
            ar & edges;
        }
    } // namespace serialization
} // namespace boost

class serializedData
{
    vector<int> serializedVector;
    int numRow;
    vector<int> numCols;
};


// 将嵌套的向量转换为一维向量进行序列化
vector<int> serialize(const vector<vector<int>>& nestedVector, int& numRow, vector<int>& numCols) {
    vector<int> serializedVector;
    numRow = nestedVector.size();
    for (const auto& innerVector : nestedVector) {
        numCols.push_back(innerVector.size());
        serializedVector.insert(serializedVector.end(), innerVector.begin(), innerVector.end());
    }
    return serializedVector;
}

// 从一维向量中反序列化恢复嵌套的向量
vector<vector<int>> deserialize(const vector<int>& serializedVector, const int& numRow, const vector<int>& numCols) {
    vector<vector<int>> nestedVector(numRow);
    size_t offset = 0;
    for (int i = 0; i < numRow; ++i) {
        size_t rowSize = numCols[i];
        nestedVector[i].insert(nestedVector[i].end(), serializedVector.begin() + offset, serializedVector.begin() + offset + rowSize);
        offset += rowSize;
    }
    return nestedVector;
}

int main(int argc, char* argv[]) {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    // MPI_Init(&argc, &argv);
    int my_rank, num_process;
    my_rank = world.rank();
    num_process = world.size();
    // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    int start_node_id, end_node_id;
    int num_vex = 0;
    double damping = 0.85;
    double tol = 1e-6;
    int max_iter = 100000;
    double start;
    unordered_map<int, int> idx2id;
    vector<vector<int>> edges;
    int num_Row;
    vector<int> num_Cols;
    int nodes_per_process;
    if(my_rank == 0) {
        ifstream infile;
        int u, v;
        unordered_map<int, int> um;
        string path = argv[1]; // "web-Google.txt"或"Data.txt"
        cout << "Processing " << path << "..." << endl;
        infile.open(path);
        int maxn = 0;
        int minn = 0x7fffffff;
        int num_edge = 0;
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
        edges = build_edges(num_vex ,um, path);
        cout<<"damping:" << damping <<"，tol:"<<tol<<endl;

        // 广播节点个数
        boost::mpi::broadcast(world, num_vex, 0);
        // MPI_Bcast(&num_vex, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // 广播edges
        vector<int> serializedData = serialize(edges, num_Row, num_Cols);

        boost::mpi::broadcast(world, serializedData, 0);
        boost::mpi::broadcast(world, num_Row, 0);
        boost::mpi::broadcast(world, num_Cols, 0);

        // 计算每个进程负责的节点范围 [start_node_id, end_node_id]
        nodes_per_process = num_vex / num_process;
        for(int i = 1; i < num_process - 1; i++) {
            start_node_id = i * nodes_per_process;
            end_node_id = start_node_id + nodes_per_process - 1;
            world.send(i, 1, start_node_id);
            world.send(i, 2, end_node_id);
            // boost::mpi::send(world, start_node_id, i, 0);
            // boost::mpi::send(world, end_node_id, i, 0);
            // MPI_Send(&start_node_id, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            // MPI_Send(&end_node_id, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
        }
        start_node_id = (num_process - 1) * nodes_per_process;
        end_node_id = num_vex - 1;
        world.send(num_process - 1, 1, start_node_id);
        world.send(num_process - 1, 2, end_node_id);
        // boost::mpi::send(world, start_node_id, num_process - 1, 0);
        // boost::mpi::send(world, end_node_id, num_process - 1, 0);
        // MPI_Send(&start_node_id, 1, MPI_INT, num_process - 1, 1, MPI_COMM_WORLD);
        // MPI_Send(&end_node_id, 1, MPI_INT, num_process - 1, 2, MPI_COMM_WORLD);
        // 主进程负责的节点范围:
        start_node_id = 0;
        end_node_id = nodes_per_process - 1;
    } else {
        // 接收节点个数
        boost::mpi::broadcast(world, num_vex, 0);
        // MPI_Bcast(&num_vex, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // 接收edges
        std::vector<int> serializedData;
        boost::mpi::broadcast(world, serializedData, 0);
        boost::mpi::broadcast(world, num_Row, 0);
        boost::mpi::broadcast(world, num_Cols, 0);
        edges = deserialize(serializedData, num_Row, num_Cols);
        // 接收每个进程负责的节点范围 [start_node_id, end_node_id]
        world.recv(0, 1, start_node_id);
        world.recv(0, 2, end_node_id);
        // boost::mpi::recv(world, 0, 0, start_node_id);
        // boost::mpi::recv(world, 0, 0, end_node_id);
        // MPI_Recv(&start_node_id, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // MPI_Recv(&end_node_id, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        nodes_per_process = end_node_id - start_node_id + 1;
    }
    // 同步所有进程
    world.barrier();
    if(my_rank == 0) {
        start = MPI_Wtime();
    }
    vector<int> out_degrees = num_Cols;
    vector<double> scores(nodes_per_process);
    // 计算 Pagerank
    int iter = 0;
    while(iter < max_iter) {
        vector<double> contribution(num_vex, 0);
        iter++;
        for(int i = start_node_id; i <= end_node_id; i++) {
            for(auto link : edges[i]) {
                contribution[link] += scores[i%nodes_per_process] / out_degrees[i];
            }
        }
        boost::mpi::all_reduce(world, boost::mpi::inplace_t<double *>(contribution.data()),  contribution.size(), std::plus<double>()); // boost.MPI原地版本
        // vector<double> newpr_reduced;
        // newpr_reduced.resize(newpr.size());
        // boost::mpi::all_reduce(world, newpr.data(), num_vex, newpr_reduced.data(), std::plus<double>()); // boost.MPI非原地版本，输出到newpr_reduced
        // MPI_Allreduce(MPI_IN_PLACE, newpr.data(), num_vex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // MPI原地版本
        // 计算新的 Pagerank 值并判断是否收敛
        double err= 0;
        for(int i = start_node_id; i <= end_node_id; i++) {
            double old_score = scores[i%nodes_per_process];
            scores[i%nodes_per_process] = (1 - damping) / num_vex + damping * contribution[i];
            err += fabs(scores[i%nodes_per_process] - old_score);
        }
        boost::mpi::all_reduce(world, boost::mpi::inplace_t<double >(err), std::plus<double>());
        if(err < tol) {
            if(my_rank == 0) {
                cout << "iter: "<< iter << ", err: "<< err << endl;
            }
            break;
        }
    }
    vector<double> gathered_scores(num_vex);
    // 计算每个进程的scores数组长度放入一个列表
    vector<int> recvcounts(num_process);
    boost::mpi::gather(world, nodes_per_process, recvcounts, 0);
    // MPI_Gather(&nodes_per_process, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // 汇总各个进程的scores数组
    boost::mpi::gatherv(world, scores.data(), nodes_per_process, gathered_scores.data(), recvcounts, 0);
    // 输出结果
    if(my_rank == 0) {
        double finish = MPI_Wtime();
        double seconds = finish - start;
        vector<pair<int, double> > indexed_scores;
        for(int i=0; i<gathered_scores.size(); ++i) {
            indexed_scores.push_back({i, gathered_scores[i]});
        }
        sort(indexed_scores.begin(), indexed_scores.end(), cmp);
        print_pagerank_scores(indexed_scores, idx2id, 10);
        cout<< "用时" << seconds << "s" << endl;
    }
    return 0;
}

/* 
//shell
mpic++ -o2  -std=c++17 -I/opt/homebrew/Cellar/boost/1.76.0/include PageRank_MPI.cpp  -L /opt/homebrew/Cellar/boost-mpi/1.76.0/lib -o mpi -lboost_mpi -lboost_serialization -L/opt/homebrew/Cellar/boost/1.76.0/lib

// MPI + small Dataset:
mpirun -np 4 ./mpi Data.txt

// MPI + large Dataset:
mpirun -np 4 ./mpi web-Google.txt
*/

/*
//debug 版本
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <mpi.h>
#include "helper_utils.h"
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/mpi.hpp>
using std::vector, std::max, std::min;

namespace boost {
    namespace serialization {
        template<class Archive>
        void serialize(Archive& ar, std::vector<std::vector<int>>& edges, const unsigned int version)
        {
            ar & edges;
        }
    } // namespace serialization
} // namespace boost

class serializedData
{
    vector<int> serializedVector;
    int numRow;
    vector<int> numCols;
};


// 将嵌套的向量转换为一维向量进行序列化
vector<int> serialize(const vector<vector<int>>& nestedVector, int& numRow, vector<int>& numCols) {
    vector<int> serializedVector;
    numRow = nestedVector.size();
    for (const auto& innerVector : nestedVector) {
        numCols.push_back(innerVector.size());
        serializedVector.insert(serializedVector.end(), innerVector.begin(), innerVector.end());
    }
    return serializedVector;
}

// 从一维向量中反序列化恢复嵌套的向量
vector<vector<int>> deserialize(const vector<int>& serializedVector, const int& numRow, const vector<int>& numCols) {
    vector<vector<int>> nestedVector(numRow);
    size_t offset = 0;
    for (int i = 0; i < numRow; ++i) {
        size_t rowSize = numCols[i];
        nestedVector[i].insert(nestedVector[i].end(), serializedVector.begin() + offset, serializedVector.begin() + offset + rowSize);
        offset += rowSize;
    }
    return nestedVector;
}

int main(int argc, char* argv[]) {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    // MPI_Init(&argc, &argv);
    int my_rank, num_process;
    my_rank = world.rank();
    num_process = world.size();
    // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    cout << "Rank: " << my_rank << "  Num_process: " << num_process << endl;
    int start_node_id, end_node_id;
    int num_vex = 0;
    double damping = 0.85;
    double tol = 1e-6;
    int max_iter = 100000;
    double start;
    unordered_map<int, int> idx2id;
    vector<vector<int>> edges;
    int num_Row;
    vector<int> num_Cols;
    int nodes_per_process;
    if(my_rank == 0) {
        cout<<"进入进程"<<my_rank<<endl;
        ifstream infile;
        int u, v;
        unordered_map<int, int> um;
        string path = "Data.txt";
        infile.open(path);
        int maxn = 0;
        int minn = 0x7fffffff;
        int num_edge = 0;
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
        edges = build_edges(num_vex ,um, path);
        cout<<"damping:" << damping <<"，tol:"<<tol<<endl;

        // TODO 调试信息
        cout<<"主进程建图成功，edges[0] from " <<my_rank<<endl;
        for(auto i: edges[0]) {
            cout<<i<<" ";
        }
        cout<<endl;
        cout<<"edges[1] from " <<my_rank<<endl;
        for(auto i: edges[1]) {
            cout<<i<<" ";
        }
        cout<<endl;
        // 广播节点个数
        boost::mpi::broadcast(world, num_vex, 0);
        // MPI_Bcast(&num_vex, 1, MPI_INT, 0, MPI_COMM_WORLD);
        cout<< my_rank << "进程broadcast num_vex成功，此时num_vex为 "<< num_vex << endl;
        // 广播edges
        vector<int> serializedData = serialize(edges, num_Row, num_Cols);
        cout<<"序列化后的1D vector的前20个元素:" <<endl;
        for(int i=0; i<20; ++i) {
            cout<<serializedData[i]<<" ";
        }
        cout << endl;
        cout << "num_Row:" << num_Row << endl;
        cout << "num_Cols的前3个元素:" << endl;
        for(int i=0; i<3;++i) {
            cout<<num_Cols[i]<<" ";
        }
        cout << endl;
        boost::mpi::broadcast(world, serializedData, 0);
        boost::mpi::broadcast(world, num_Row, 0);
        boost::mpi::broadcast(world, num_Cols, 0);
        cout<< my_rank << "进程broadcast edges成功，此时edges[0]为:"<< endl;
        for(auto i: edges[0]) {
            cout<<i<<" ";
        }
        cout << endl;
        // 计算每个进程负责的节点范围 [start_node_id, end_node_id]
        nodes_per_process = num_vex / num_process;
        for(int i = 1; i < num_process - 1; i++) {
            start_node_id = i * nodes_per_process;
            end_node_id = start_node_id + nodes_per_process - 1;
            world.send(i, 1, start_node_id);
            world.send(i, 2, end_node_id);
            // boost::mpi::send(world, start_node_id, i, 0);
            // boost::mpi::send(world, end_node_id, i, 0);
            // MPI_Send(&start_node_id, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            // MPI_Send(&end_node_id, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
        }
        start_node_id = (num_process - 1) * nodes_per_process;
        end_node_id = num_vex - 1;
        world.send(num_process - 1, 1, start_node_id);
        world.send(num_process - 1, 2, end_node_id);
        // boost::mpi::send(world, start_node_id, num_process - 1, 0);
        // boost::mpi::send(world, end_node_id, num_process - 1, 0);
        // MPI_Send(&start_node_id, 1, MPI_INT, num_process - 1, 1, MPI_COMM_WORLD);
        // MPI_Send(&end_node_id, 1, MPI_INT, num_process - 1, 2, MPI_COMM_WORLD);
        // 主进程负责的节点范围:
        cout<< my_rank << "进程send to "<< num_process - 1 << "进程 start_node_id、end_node_id成功，此时start_node_id、end_node_id为:"<< start_node_id << " " << end_node_id << endl;
        start_node_id = 0;
        end_node_id = nodes_per_process - 1;
        cout<< my_rank << "进程负责的start_node_id、end_node_id为: "<< start_node_id << " " << end_node_id << endl;
        start = MPI_Wtime();
    } else {
        cout<<"进入进程"<<my_rank<<endl;
        // 接收节点个数
        boost::mpi::broadcast(world, num_vex, 0);
        // MPI_Bcast(&num_vex, 1, MPI_INT, 0, MPI_COMM_WORLD);
        cout<< my_rank << "进程broadcast num_vex成功，此时num_vex为 "<< num_vex<<endl;
        // 接收edges
        std::vector<int> serializedData;
        boost::mpi::broadcast(world, serializedData, 0);
        boost::mpi::broadcast(world, num_Row, 0);
        boost::mpi::broadcast(world, num_Cols, 0);
        edges = deserialize(serializedData, num_Row, num_Cols);
        cout<< my_rank << "进程broadcast edges成功，此时edges[0]为:"<< endl;
        for(auto i: edges[0]) {
            cout<<i<<" ";
        }
        cout << endl;
        cout << "num_Row:" << num_Row << endl;
        cout << "num_Cols的前3个元素:" << endl;
        for(int i=0; i<3;++i) {
            cout<<num_Cols[i]<<" ";
        }
        cout<<endl;
        // 接收每个进程负责的节点范围 [start_node_id, end_node_id]
        world.recv(0, 1, start_node_id);
        world.recv(0, 2, end_node_id);
        // boost::mpi::recv(world, 0, 0, start_node_id);
        // boost::mpi::recv(world, 0, 0, end_node_id);
        // MPI_Recv(&start_node_id, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // MPI_Recv(&end_node_id, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        nodes_per_process = end_node_id - start_node_id + 1;
        cout<< my_rank << "进程recv start_node_id、end_node_id成功，此时start_node_id、end_node_id为:"<< start_node_id << " " << end_node_id << endl;
    }
    vector<int> out_degrees = num_Cols;
    vector<double> scores(nodes_per_process);
    // 计算 Pagerank
    int iter = 0;
    while(iter < max_iter) {
        vector<double> contribution(num_vex, 0);
        iter++;
        for(int i = start_node_id; i <= end_node_id; i++) {
            for(auto link : edges[i]) {
                contribution[link] += scores[i%nodes_per_process] / out_degrees[i];
            }
        }
        cout<< my_rank << "进程all_reduce之前，此时contribution的前20个元素为:" << endl;
        for(int i=0; i<10; ++i) {
            cout<<contribution[i]<<" ";
        }
        cout << endl;
        boost::mpi::all_reduce(world, boost::mpi::inplace_t<double *>(contribution.data()),  contribution.size(), std::plus<double>()); // boost.MPI原地版本
        // vector<double> newpr_reduced;
        // newpr_reduced.resize(newpr.size());
        // boost::mpi::all_reduce(world, newpr.data(), num_vex, newpr_reduced.data(), std::plus<double>()); // boost.MPI非原地版本，输出到newpr_reduced
        // MPI_Allreduce(MPI_IN_PLACE, newpr.data(), num_vex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // MPI原地版本
        if(my_rank == 0) {
            cout<< my_rank << "进程all_reduce成功，此时contribution的前10个元素为:" << endl;
            for(int i=0; i<10; ++i) {
                cout<<contribution[i]<<" "; 
            }
            cout << endl;
            cout<< my_rank << "进程all_gather前scores: " << endl;
            for(int i=start_node_id; i<=start_node_id+10; ++i) {
                cout<<scores[i%nodes_per_process]<<" "; 
            }
            cout << endl;
        }
        // 计算新的 Pagerank 值并判断是否收敛
        double err= 0;
        for(int i = start_node_id; i <= end_node_id; i++) {
            double old_score = scores[i%nodes_per_process];
            scores[i%nodes_per_process] = (1 - damping) / num_vex + damping * contribution[i];
            err += fabs(scores[i%nodes_per_process] - old_score);
        }
        cout<< my_rank << "进程all_reduce前err: " << err << endl;
        boost::mpi::all_reduce(world, boost::mpi::inplace_t<double >(err), std::plus<double>());
        if(my_rank == 0) {
            cout<< my_rank << "进程all_reduce后err: " << err << endl;
        }
        if(err < tol) {
            if(my_rank == 0) {
                cout << "iter: "<< iter << ", err: "<< err << endl;
            }
            break;
        }
    }
    vector<double> gathered_scores(num_vex);
    cout<< my_rank << "进程nodes_per_process:" << nodes_per_process << endl;
    // 计算每个进程的scores数组长度放入一个列表
    vector<int> recvcounts(num_process);
    boost::mpi::gather(world, nodes_per_process, recvcounts, 0);
    // MPI_Gather(&nodes_per_process, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(my_rank == 0) {
        cout<< my_rank << "进程汇总各进程scores数组长度的列表:" << endl;
        for(auto i:recvcounts) {
            cout << i << " ";
        }
        cout<<endl;
        //汇总后总长度其实就是num_vex
    }
    // 同步所有进程
    world.barrier();
    // MPI_Barrier(MPI_COMM_WORLD);
    // 汇总各个进程的scores数组
    boost::mpi::gatherv(world, scores.data(), nodes_per_process, gathered_scores.data(), recvcounts, 0);
    // 输出结果
    if(my_rank == 0) {
        double finish = MPI_Wtime();
        double seconds = finish - start;
        vector<pair<int, double> > indexed_scores;
        for(int i=0; i<gathered_scores.size(); ++i) {
            indexed_scores.push_back({i, gathered_scores[i]});
        }
        sort(indexed_scores.begin(), indexed_scores.end(), cmp);
        print_pagerank_scores(indexed_scores, idx2id, 100);
        cout<< "用时" << seconds << "s" << endl;
    }
    // MPI_Finalize();
    return 0;
}
*/



