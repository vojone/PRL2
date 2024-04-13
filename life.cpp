/**
 * @file life.cpp
 * @author Vojtech Dvorak (xdvora3o@fit.vutbr.cz)
 *
 */

//#include <mpi.h>

#include <iostream>
#include <fstream>
#include <vector>



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
            std::cerr << c << std::endl;
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
    if(argc < 3) {
        std::cerr << "USAGE: ./life <path-to-board> <it-num> " << std::endl;
        return 2;
    }

    const char *board_path = argv[1];
    size_t iteration_number = std::stoi(argv[2]);

    int rank = 0, size;
    // MPI_Init (&argc, &argv);

    // // Get rank and total number of processors
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank == 0) { // The first (input) processor

        std::vector<std::vector<uint8_t>> board;

        int ret = parse_board(argv[1], board);
        if(ret) {
            return ret;
        }
    }

    return 0;
}
