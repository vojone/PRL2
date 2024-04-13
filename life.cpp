/**
 * @file life.cpp
 * @author Vojtech Dvorak (xdvora3o@fit.vutbr.cz)
 *
 */

#include <mpi.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define DEBUG


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




int main(int argc, char **argv) {
    MPI_Init (&argc, &argv);

    int rank = 0, size, ret = 0;

    // Get rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm comm;
    int dims[1] = {size}, periodic[1] = {1}, coords[1];
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periodic, 0, &comm);


    std::vector<uint8_t> board, chunk;
    size_t row_len = 0, rows_per_process = 0, row_num = 0, remaining_rows_num = 0;
    if(rank == 0) { // The root processor processor
        if(argc < 3) {
            std::cerr << "USAGE: ./life <path-to-board> <it-num> " << std::endl;
            return 2;
        }

        const char *board_path = argv[1];
        size_t iteration_number = std::stoi(argv[2]);


        ret = parse_board(argv[1], board, row_len, row_num);
        if(ret) {
            MPI_Finalize();
            return ret;
        }

        rows_per_process = (size_t)ceil(row_num / ((float)size));
        remaining_rows_num = row_num - rows_per_process * (size - 1);
        LOG(rank, "Board size: %dx%d", row_len, row_num);
        LOG(rank, "Rows per process: %d", rows_per_process);
        LOG(rank, "Remaining rows: %d", remaining_rows_num);
    }
    
    MPI_Bcast(&row_len, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rows_per_process, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    size_t remaining_rows_num = row_num - rows_per_process * (size - 1);
    
    int displs[size], send_counts[size];
    send_counts[0] = 0;
    displs[0] = 0;
    for(size_t i = 1; i < size; ++i) {
        send_counts[i] = row_len * rows_per_process;
        displs[i] = row_len * remaining_rows_num  + (i - 1) * row_len * rows_per_process;
    }

    if(rank != 0) {
        chunk.resize(row_len * rows_per_process);
    }
    else {
        chunk.assign(board.begin(), board.begin() + row_len * remaining_rows_num);
    }
    
    MPI_Scatterv(board.data(), send_counts, displs, MPI_CHAR, chunk.data(), row_len * rows_per_process, MPI_CHAR, 0, MPI_COMM_WORLD);



    MPI_Finalize();

    return ret;
}
