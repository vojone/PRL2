/**
 * @file life.cpp
 * @author Vojtech Dvorak (xdvora3o@fit.vutbr.cz)
 *
 */

#include <mpi.h>

#include <iostream>
#include <fstream>
#include <vector>

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


int parse_board(const char* input_file_path, std::vector<std::vector<uint8_t>> &board) {
    std::ifstream input_file(input_file_path, std::ios::in);

    board.clear();
    board.push_back(std::vector<uint8_t>());

    char c;
    size_t line_i = 0;
    while(!input_file.eof()) {
        input_file.get(c);
        if(input_file.eof()) {
            break;
        }
    
        if(c == '\n' || c == '\r') {
            if(c == '\r') {
                input_file.get(c);
            }

            if(board[line_i].size() != board.front().size()) {
                std::cerr << "line " << (line_i + 1);
                std::cerr << ": The input board must have a rectangular shape!" << std::endl;
                return 1;
            }

            line_i++;
            board.push_back(std::vector<uint8_t>({}));
        }
        else {
            if(c == '1') {
                board[line_i].push_back(1);
            }
            else {
                board[line_i].push_back(0);
            }
        }
    }

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


    std::vector<std::vector<uint8_t>> my_cells;
    if(rank == 0) { // The first (input) processor
        if(argc < 3) {
            std::cerr << "USAGE: ./life <path-to-board> <it-num> " << std::endl;
            return 2;
        }

        const char *board_path = argv[1];
        size_t iteration_number = std::stoi(argv[2]);

        std::vector<std::vector<uint8_t>> board;

        ret = parse_board(argv[1], board);

    }
    else {

    }

    int left, right;
    MPI_Cart_shift(comm, 0, 1, &left, &right);
    LOG(rank, "L: %d R: %d", left, right);

    MPI_Finalize();

    return ret;
}
