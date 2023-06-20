// // #include <stdio.h>
// // #include <mpi.h>
// // using namespace std;
// // int main(int argc, char* argv[]) {
// //     int rank, size;

// //     // 初始化 MPI
// //     MPI_Init(&argc, &argv);

// //     // 获取当前进程的Rank
// //     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// //     // 获取进程总数
// //     MPI_Comm_size(MPI_COMM_WORLD, &size);

// //     // 打印 Hello World 消息
// //     printf("Hello World from process %d of %d\n", rank, size);

// //     // 终止 MPI
// //     MPI_Finalize();

// //     return 0;
// // }

// #include <stdio.h>
// #include <string.h>
// #include "mpi.h"
// int main(int argc, char* argv[]) {
//     int numprocs, myid, source;
//     MPI_Status status;
//     char message[100];
//     MPI_Init(&argc, &argv);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//     MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//     if (myid != 0) {  //非0号进程发送消息
//         strcpy(message, "Hello World!");
//         MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 99,
//             MPI_COMM_WORLD);
//     }
//     else {   // myid == 0，即0号进程接收消息
//         for (source = 1; source < numprocs; source++) {
//             MPI_Recv(message, 100, MPI_CHAR, source, 99,
//                 MPI_COMM_WORLD, &status);
//             printf("接收到第%d号进程发送的消息：%s\n", source, message);
//         }
//     }
//     MPI_Finalize();
//     return 0;
// } /* end main */




// namespace boost {
//     namespace serialization {
//         template<class Archive>
//         void serialize(Archive& ar, std::vector<std::vector<int>>& edges, const unsigned int version)
//         {
//             ar & edges;
//         }
//     } // namespace serialization
// } // namespace boost

// int main(int argc, char* argv[]) {
//     // boost::mpi::environment env(argc, argv);
//     boost::mpi::communicator world;
//     MPI_Init(&argc, &argv);
//     int my_rank, num_process;
//     my_rank = world.rank();
//     num_process = world.size();
//     // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//     // MPI_Comm_size(MPI_COMM_WORLD, &num_process);
//     cout << "Rank: " << my_rank << "  Num_process: " << num_process << endl;
//     int start_node_id, end_node_id;
//     int num_vex = 0;
//     double damping = 0.85;
//     double tol = 1e-6;
//     int max_iter = 1000;
//     double start;
//     unordered_map<int, int> idx2id;
//     vector<vector<int>> edges;
//     int nodes_per_process;
//     if(my_rank == 0) {
//         cout<<"haha";
//         ifstream infile;
//         int u, v;
//         unordered_map<int, int> um;
//         string path = "Data.txt";
//         infile.open(path);
//         int maxn = 0;
//         int minn = 0x7fffffff;
//         int num_edge = 0;
//         int cnt = 0;
//         while(infile>>u>>v) {
//             if(um.find(u)==um.end()) {
//                 um[u] = num_vex++;
//             }
//             if(um.find(v)==um.end()) {
//                 um[v] = num_vex++;
//             }
//             ++num_edge;
//             maxn = max(max(maxn, u),v);
//             minn = min(min(minn, u),v);
//         }
//         for(auto kv : um) {
//             idx2id[kv.second] = kv.first;
//         }
//         ofstream outfile("idx2id.txt");
//         if (!outfile) {
//             cout << "Error opening file." << endl;
//             exit(-1);
//         }
//         for (auto it = idx2id.begin(); it != idx2id.end(); it++) {
//             outfile << it->first << " " << it->second << endl;
//         }
//         outfile.close();
//         cout<<"最大节点为："<<maxn<<endl;
//         cout<<"最小节点为："<<minn<<endl;
//         cout<<"节点个数："<<num_vex<<"，边个数："<<num_edge<<endl;
//         infile.close();
//         edges = build_edges(num_vex ,um, path);
//         cout<<"damping:" << damping <<"，tol:"<<tol<<endl;
//         // 广播节点个数
//         boost::mpi::broadcast(world, num_vex, 0);
//         MPI_Send(&num_vex, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
//         cout<<'23123';
//         MPI_Bcast(&num_vex, 1, MPI_INT, 0, MPI_COMM_WORLD);
//         cout<<12312;
//         // 广播edges
//         boost::mpi::broadcast(world, edges, 0);
//         // 计算每个进程负责的节点范围 [start_node_id, end_node_id]
//         cout<<12312;
//         nodes_per_process = num_vex / num_process;
//         for(int i = 1; i < num_process - 1; i++) {
//             start_node_id = i * nodes_per_process;
//             end_node_id = start_node_id + nodes_per_process - 1;
//             world.send(i, 1, start_node_id);
//             world.send(i, 2, end_node_id);
//             // boost::mpi::send(world, start_node_id, i, 0);
//             // boost::mpi::send(world, end_node_id, i, 0);
//             // MPI_Send(&start_node_id, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
//             // MPI_Send(&end_node_id, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
//         }
//         start_node_id = (num_process - 1) * nodes_per_process;
//         end_node_id = num_vex - 1;
//         world.send(num_process - 1, 1, start_node_id);
//         world.send(num_process - 1, 2, end_node_id);
//         // boost::mpi::send(world, start_node_id, num_process - 1, 0);
//         // boost::mpi::send(world, end_node_id, num_process - 1, 0);
//         // MPI_Send(&start_node_id, 1, MPI_INT, num_process - 1, 1, MPI_COMM_WORLD);
//         // MPI_Send(&end_node_id, 1, MPI_INT, num_process - 1, 2, MPI_COMM_WORLD);
//         // 主进程负责的节点范围:
//         start_node_id = 0;
//         end_node_id = nodes_per_process - 1;
//         start = MPI_Wtime();
//     } else {
//         cout<<"nihao";
//         MPI_Recv(&num_vex, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         cout<<"nihao2";
//         // 接收节点个数
//         // boost::mpi::broadcast(world, num_vex, 0);
//         MPI_Bcast(&num_vex, 1, MPI_INT, 0, MPI_COMM_WORLD);
//         cout<<'asd';
//         // 接收edges
//         boost::mpi::broadcast(world, edges, 0);
//         // 接收每个进程负责的节点范围 [start_node_id, end_node_id]
//         world.recv(0, 1, start_node_id);
//         world.recv(0, 2, end_node_id);
//         // boost::mpi::recv(world, 0, 0, start_node_id);
//         // boost::mpi::recv(world, 0, 0, end_node_id);
//         // MPI_Recv(&start_node_id, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         // MPI_Recv(&end_node_id, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         nodes_per_process = end_node_id - start_node_id + 1;
//         // 接收每个节点指向的节点集合
//         vector<int> send_links;
//         vector<int> recv_links(num_vex * num_vex);
//         for(int i = start_node_id; i <= end_node_id; i++) {
//             out_degrees[i] = edges[i].size();
//             for(auto it : edges[i]) {
//                 send_links.push_back(it);
//             }
//         }
//         boost::mpi::all_gather(world, out_degrees.data() + start_node_id, nodes_per_process, out_degrees.data());
//         boost::mpi::all_gather(world, send_links.data(), send_links.size(), recv_links.data());
//         // MPI_Allgather(out_degrees.data() + start_node_id, nodes_per_process, MPI_INT, out_degrees.data(), nodes_per_process, MPI_INT, MPI_COMM_WORLD);
//         // MPI_Allgather(send_links.data(), send_links.size(), MPI_INT, recv_links.data(), num_vex * nodes_per_process, MPI_INT, MPI_COMM_WORLD);
//         vector<vector<int>> recv_links_list(num_vex);
//         for (int i = 0; i < num_vex; i++) {
//             for (int j = 0; j < out_degrees[i]; j++) {
//                 recv_links_list[i].push_back(recv_links[i * num_vex + j]);
//             }
//         }
//     }
//     vector<int> out_degrees(num_vex);
//     for(int i = start_node_id; i <= end_node_id; i++) {
//         out_degrees[i] = edges[i].size();
//     }
//     boost::mpi::all_gather(world, out_degrees.data() + start_node_id, nodes_per_process, out_degrees.data());
//     // 初始化 Pagerank 值和每个节点的出链数量
//     vector<double> scores(nodes_per_process, 1.0 / num_vex);
//     vector<double> newpr(num_vex, (1.0 - damping) / num_vex);
//     // 计算 Pagerank
//     int iteration = 0;
//     while(iteration < max_iter) {
//         iteration++;
//         // 计算每个节点的贡献值
//         for(int i = start_node_id; i <= end_node_id; i++) {
//             for(auto link : edges[i]) {
//                 newpr[link] += damping * (scores[i] / out_degrees[i]);
//             }
//         }
//         boost::mpi::all_reduce(world, newpr.data(), num_vex, newpr.data(), std::plus<int>());
//         // MPI_Allreduce(MPI_IN_PLACE, contributions.data(), num_vex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//         // 计算新的 Pagerank 值并判断是否收敛
//         double local_err = 0;
//         double err;
//         for(int i = start_node_id; i <= end_node_id; i++) {
//             local_err += fabs(scores[i] - newpr[i]);
// 			scores[i] = newpr[i];
// 			newpr[i] = (1.0 - damping) / num_vex;
//                 // double old_pagerank = scores[i];
//                 // scores[i] = (1 - damping) / num_vex + damping * contributions[i];
//                 // sum += scores[i];
//                 // if(fabs(scores[i] - old_pagerank) < tol) {
//                 //     is_converged = true;
//                 // }
//         }
//         boost::mpi::all_reduce(world, local_err, err, std::plus<int>());
//         if(err < tol) {
//             if(my_rank == 0) {
//                 cout<<"err:"<<err<<endl;
//             }
//             break;
//         }
//         // MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         // for(int i = start_node_id; i <= end_node_id; i++) {
//         //     scores[i] += (1 - global_sum) / num_vex;
//         // }
//         // 同步 Pagerank 值
//         MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
//                         scores.data(), num_vex / num_process, MPI_DOUBLE, MPI_COMM_WORLD);
//     }

//     // 输出结果
//     if(my_rank == 0) {
//         double finish = MPI_Wtime();
//         double seconds = finish - start;
//         vector<pair<int, double> > indexed_scores;
//         for(int i=0; i<scores.size(); ++i) {
//             indexed_scores.push_back({i, scores[i]});
//         }
//         sort(indexed_scores.begin(), indexed_scores.end(), cmp);
//         print_pagerank_scores(indexed_scores, idx2id, 100);
//         cout<< "用时" << seconds << "s" << endl;
//     }
//     MPI_Finalize();
//     return 0;
// }

#include <iostream>
#include <vector>
#include <boost/mpi.hpp>

int main() {
    boost::mpi::environment env;
    boost::mpi::communicator world;

    // int my_vector[5] = {1, 2, 3, 4, 5};
    // int result_vector[5];
    std::vector<int> my_vector = {1, 2, 3, 4, 5};
    std::vector<int> result_vector;
    result_vector.resize(my_vector.size());
    boost::mpi::all_reduce(world, my_vector.data(), 5, result_vector.data(), std::plus<int>());

    if (world.rank() == 0) {
        std::cout << "Result vector: ";
        for (int num : result_vector) {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}






