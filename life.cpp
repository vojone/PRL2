/**
 * @file life.cpp
 * @author Vojtech Dvorak (xdvora3o@fit.vutbr.cz)
 *
 */

#include <mpi.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>

#define DEBUG


#define HISTORY_LEN 2

typedef enum { BOT, TOP, NEIGH_NUM } neighbour_index_t;


#ifdef DEBUG
    #define LOG(rank, ...)\
        fprintf(stderr, "Rank: %d  ", rank);\
        fprintf(stderr, __VA_ARGS__);\
        fprintf(stderr, "\n");\
        fflush(stderr);
#else
    #define LOG(...) ""
#endif


int parse_board(const char* input_file_path, std::vector<uint8_t> &board, size_t &row_len, size_t &row_num) {
    std::ifstream input_file(input_file_path, std::ios::in);

    board.clear();
    row_len = -1;

    char c;
    size_t line_i = 0, cur_row_len = 0;
    while(!input_file.eof()) {
        input_file.get(c);
        if(input_file.eof()) {
            break;
        }
    
        if(c == '\n' || c == '\r') {
            if(c == '\r') {
                input_file.get(c);
            }

            if(line_i == 0) {
                row_len = cur_row_len;
            }
            else {
                if(row_len != cur_row_len) {
                    std::cerr << "line " << (line_i + 1);
                    std::cerr << ": The input board must have a rectangular shape!" << std::endl;
                    return 1;
                }
            }

            line_i++;
            cur_row_len = 0;
        }
        else {
            if(c == '1') {
                board.push_back(1);
            }
            else {
                board.push_back(0);
            }

            cur_row_len++;
        }
    }

    row_num = line_i;
    return 0;
}


void print_board(int rank, std::vector<uint8_t> board, size_t row_num, size_t row_len, size_t rows_per_proc, size_t remaining_rows_num, bool print_prefixes) {
    for(size_t i = 0; i < row_num; ++i) {
        if(print_prefixes) {
            if(i < remaining_rows_num) {
                std::cout << rank << ": ";
            }
            else {
                std::cout << (i - remaining_rows_num) / rows_per_proc + 1 << ": ";
            }
        }

        for(size_t j = 0; j < row_len; ++j) {
            std::cout << (int)board[i * row_len + j];
        }

        std::cout << std::endl;
    }
}


uint8_t get_state_wa(std::vector<uint8_t> &chunk, std::vector<uint8_t> *neigh_rows, int row_len, int r, int c) {
    size_t row_num = (chunk.size() / row_len);
    if(c >= row_len) {
        c = 0;
    }
    else if(c < 0) {
        c += row_len;
    }

    if(r >= (int)row_num) {
        return neigh_rows[TOP][c + row_len];
    }
    else if(r < 0) {
        return neigh_rows[BOT][c];
    }
    else {
        return chunk[r * row_len + c];
    }
}

int main(int argc, char **argv) {
    MPI_Init (&argc, &argv);

    int rank = 0, size, ret = 0;

    // Get rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm comm;
    int dims[1] = {size}, periodic[1] = {1}, coords[1];
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periodic, 0, &comm);


    std::vector<uint8_t> board, neigh_rows[NEIGH_NUM], chunk[HISTORY_LEN];
    size_t row_len = 0, rows_per_process = 0, row_num = 0, remaining_rows_num = 0;
    size_t iteration_number;
    if(rank == 0) { // The root processor processor
        if(argc < 3) {
            std::cerr << "USAGE: ./life <path-to-board> <it-num> " << std::endl;
            return 2;
        }

        const char *board_path = argv[1];
        iteration_number = std::stoi(argv[2]);


        ret = parse_board(argv[1], board, row_len, row_num);
        if(ret) {
            MPI_Finalize();
            return ret;
        }

        rows_per_process = (size_t)ceil(row_num / ((float)size));
        remaining_rows_num = row_num - rows_per_process * (size - 1);
        LOG(rank, "Board size: %ldx%ld", row_len, row_num);
        LOG(rank, "Rows per process: %ld", rows_per_process);
        LOG(rank, "Remaining rows: %ld", remaining_rows_num);
    }
    
    MPI_Bcast(&row_len, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rows_per_process, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iteration_number, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    
    int displs[size], send_counts[size];
    send_counts[0] = 0;
    displs[0] = 0;
    for(size_t i = 1; i < size; ++i) {
        send_counts[i] = row_len * rows_per_process;
        displs[i] = row_len * remaining_rows_num  + (i - 1) * row_len * rows_per_process;
    }

    if(rank != 0) {
        chunk[0].resize(row_len * rows_per_process);
    }
    else {
        chunk[0].assign(board.begin(), board.begin() + row_len * remaining_rows_num);
    }
    
    MPI_Scatterv(board.data(), send_counts, displs, MPI_CHAR, chunk[0].data(), row_len * rows_per_process, MPI_CHAR, 0, MPI_COMM_WORLD);
    chunk[1].assign(chunk[0].begin(), chunk[0].end());


    neigh_rows[BOT].resize(NEIGH_NUM * row_len);
    neigh_rows[TOP].resize(NEIGH_NUM * row_len);

    int cur_state = 0, next_state = 1;
    for(size_t i = 0; i < iteration_number; ++i) {
        uint8_t* last_row_ptr = &(chunk[cur_state].data()[chunk[cur_state].size() - row_len]);

        MPI_Neighbor_allgather(last_row_ptr, row_len, MPI_CHAR, neigh_rows[BOT].data(), row_len, MPI_CHAR, comm);
        MPI_Neighbor_allgather(chunk[cur_state].data(), row_len, MPI_CHAR, neigh_rows[TOP].data(), row_len, MPI_CHAR, comm);

        for(int r = 0; r < chunk[cur_state].size() / row_len; ++r) {
            for(int c = 0; c < row_len; ++c) {
                int i = r * row_len + c;

                if(rank == 2) {
                    std::cerr << "Cell " << r << ", " << c << std::endl;

                    uint8_t s;
                    s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r - 1, c - 1);
                    std::cerr << (int)s;
                    s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r - 1, c);
                    std::cerr << (int)s;
                    s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r - 1, c + 1);
                    std::cerr << (int)s;

                     std::cerr << std::endl;

                    s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r, c - 1);
                    std::cerr << (int)s;
                    // s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r, c);
                    // std::cerr << (int)s;
                    std::cerr << " ";
                    s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r, c + 1);
                    std::cerr << (int)s;

                     std::cerr << std::endl;

                    s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r + 1, c - 1);
                    std::cerr << (int)s;
                    s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r + 1, c);
                    std::cerr << (int)s;
                    s = get_state_wa(chunk[cur_state], neigh_rows, row_len, r + 1, c + 1);
                    std::cerr << (int)s;

                    std::cerr << std::endl;
                }



                chunk[next_state][i] = 1;
                // if(rank == 0) {
                
                //     std::cout << get_state_wa(chunk, neigh_rows[])
                // }


            }
        }

        if(rank == 2) {
            int left, right;
            MPI_Cart_shift(comm, 0, 1, &left, &right);
            LOG(rank, "L: %d R: %d", left, right);
            LOG(rank, "");
            for(auto &c: neigh_rows[BOT]) {
                std::cerr << (int)c;
            }
            std::cerr << std::endl;
            for(auto &c: neigh_rows[TOP]) {
                std::cerr << (int)c;
            }
            std::cerr << std::endl;
            fflush(stderr);
        }
    }


    MPI_Gatherv(chunk[0].data(), row_len * rows_per_process, MPI_CHAR, board.data(), send_counts, displs, MPI_CHAR, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        print_board(rank, board, row_num, row_len, rows_per_process, remaining_rows_num, true);
    }

    MPI_Finalize();

    return ret;
}
