/**
 * @file life.cpp
 * @author Vojtech Dvorak (xdvora3o@fit.vutbr.cz)
 *
 */

#include <mpi.h>

#include <iostream>
#include <fstream>
#include <vector>



int parse_board(const char* input_file_path, std::vector<std::vector<uint8_t>> &board) {
    std::ifstream input_file(input_file_path, std::ios::in);

    board.clear();
    board.insert(board.begin(), std::vector<uint8_t>());

    char c;
    size_t line_cnt = 1;
    while(!input_file.eof()) {
        input_file.get(c);
        if(input_file.eof()) {
            break;
        }
        else if(c == '\n' || c == '\r') {
            if(c == '\r') {
                input_file.get(c);
            }

            if(board.size() > 1 && board[board.size() - 1].size() != board[board.size() - 2].size()) {
                std::cerr << "line " << line_cnt << ": The input board must have rectangular shape!" << std::endl;
                return 1;
            }

            line_cnt++;
            board.insert(board.begin(), std::vector<uint8_t>());
        }
        else {
            if(c == '1') {
                board[board.size() - 1].insert(board[board.size() - 1].begin(), 1);
            }
            else {
                board[board.size() - 1].insert(board[board.size() - 1].begin(), 0);
            }
        }
    }

    return 0;
}




int main(int argc, char **argv) {
    if(argc < 2) {
        std::cerr << "USAGE: ./life <path-to-board> <it-num> " << std::endl;
        return 2;
    }

    std::vector<std::vector<uint8_t>> board;

    int ret = parse_board(argv[0], board);
    if(ret) {
        return ret;
    }

    return 0;
}
