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
#define BOARD_DIM 2

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


typedef struct gol_conf {
    size_t row_num;
    size_t row_len;
    size_t rows_per_proc;
    size_t master_rows_num;
    size_t iteration_n;
} gol_conf_t;


typedef struct test {
    size_t row_num;
    size_t row_len;
    size_t rows_per_proc;
    size_t master_rows_num;
    size_t iteration_n;
} test_t;

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
            board.push_back(c == '1' ? 1 : 0);
            cur_row_len++;
        }
    }

    row_num = line_i;
    return 0;
}


void print_board(int rank, int *proc_rows, std::vector<uint8_t> board, gol_conf_t &conf, bool print_prefixes) {
    size_t p = 0;
    for(size_t i = 0, proc_row_cnt = 0; i < conf.row_num; ++i, ++proc_row_cnt) {
        if(print_prefixes) {
            while(proc_rows[p] <= proc_row_cnt) {
                proc_row_cnt = 0;
                ++p;
            }

            std::cout << p << ": ";
        }

        for(size_t j = 0; j < conf.row_len; ++j) {
            std::cout << (int)board[i * conf.row_len + j];
        }

        std::cout << std::endl;
    }
}


uint8_t get_state_wa(
    const std::vector<uint8_t> &chunk,
    const std::vector<uint8_t> *neigh_rows,
    int row_len,
    int row,
    int col)
{
    int row_num = chunk.size() / row_len;
    if(col >= row_len) {
        col = 0;
    }
    else if(col < 0) {
        col += row_len;
    }

    const uint8_t *source = chunk.data();
    if(row >= row_num) {
        row = 0;
        source = &(neigh_rows[TOP][row_len]);
    }
    else if(row < 0) {
        row = 0;
        source = &(neigh_rows[BOT][0]);
    }

    return source[row * row_len + col];
}


uint8_t get_state_closed(
    int rank,
    int size,
    const std::vector<uint8_t> &chunk,
    const std::vector<uint8_t> *neigh_rows,
    int row_len,
    int row,
    int col)
{
    int row_num = chunk.size() / row_len;
    if(col >= row_len || col < 0) {
        return 0;
    }

    const uint8_t *source = chunk.data();
    if(row >= row_num) {
        if(rank >= size - 1) {
            return 0;
        }

        row = 0;
        source = &(neigh_rows[TOP][row_len]);
    }
    else if(row < 0) {
        if(rank == 0) {
            return 0;
        }

        row = 0;
        source = &(neigh_rows[BOT][0]);
    }

    return source[row * row_len + col];
}



int main(int argc, char **argv) {
    MPI_Init (&argc, &argv);

    int rank = 0, size;

    // Get rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    std::vector<uint8_t> board, neigh_rows[NEIGH_NUM], chunk[HISTORY_LEN];
    int proc_rows[size];
    gol_conf_t conf = { .row_num = 0, .row_len = 0, .rows_per_proc = 0, .master_rows_num = 0, .iteration_n = 0 };
    if(rank == 0) { // The root processor processor
        if(argc < 3) {
            std::cerr << "USAGE: ./life <path-to-board> <it-num> " << std::endl;
            return 2;
        }

        const char *board_path = argv[1];
        conf.iteration_n = std::stoi(argv[2]);


        int ret = parse_board(argv[1], board, conf.row_len, conf.row_num);
        if(ret) {
            MPI_Finalize();
            return ret;
        }

        conf.rows_per_proc = (size_t)(conf.row_num / ((float)(size - 1)));
        size_t remaining_rows = conf.row_num; 
        for(size_t i = 0; remaining_rows >= conf.rows_per_proc && i < size - 1; ++i) {
            proc_rows[i] = conf.rows_per_proc;
            remaining_rows -= conf.rows_per_proc;
        }

        proc_rows[size - 1] = remaining_rows;

        LOG(rank, "Board size: %ldx%ld", conf.row_len, conf.row_num);
        LOG(rank, "Rows per process: %ld", conf.rows_per_proc);
        LOG(rank, "Remaining rows: %ld", proc_rows[size - 1]);
    }

    MPI_Bcast(proc_rows, size, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&conf, sizeof(conf), MPI_BYTE, 0, MPI_COMM_WORLD);

    int displs[size], proc_cells[size];
    for(size_t i = 0; i < size; ++i) {
        proc_cells[i] = conf.row_len * proc_rows[i];
        displs[i] = i * conf.row_len * conf.rows_per_proc;
    }

    chunk[0].resize(proc_rows[rank] * conf.row_len);
    
    MPI_Scatterv(board.data(), proc_cells, displs, MPI_CHAR, chunk[0].data(), proc_cells[rank], MPI_CHAR, 0, MPI_COMM_WORLD);
    chunk[1].assign(chunk[0].begin(), chunk[0].end());



    MPI_Comm chunk_holder_comm, grid;
    MPI_Comm_split(MPI_COMM_WORLD, proc_rows[rank] > 0, rank, &chunk_holder_comm);

    int grid_size;
    MPI_Comm_size(chunk_holder_comm, &grid_size);

    int dims[BOARD_DIM] = {grid_size}, periodic[1] = {1};
    MPI_Cart_create(chunk_holder_comm, 1, dims, periodic, 0, &grid);

    neigh_rows[BOT].resize(NEIGH_NUM * conf.row_len);
    neigh_rows[TOP].resize(NEIGH_NUM * conf.row_len);

    int cur_state = 0, next_state = 1;
    if(chunk[0].size()) {
        for(size_t i = 0; i < conf.iteration_n; ++i) {
            uint8_t* last_row_ptr = &(chunk[cur_state].data()[chunk[cur_state].size() - conf.row_len]);

            MPI_Neighbor_allgather(last_row_ptr, conf.row_len, MPI_CHAR, neigh_rows[BOT].data(), conf.row_len, MPI_CHAR, grid);
            MPI_Neighbor_allgather(chunk[cur_state].data(), conf.row_len, MPI_CHAR, neigh_rows[TOP].data(), conf.row_len, MPI_CHAR, grid);

            memset(chunk[next_state].data(), 0, chunk[next_state].size());
            for(int r = 0; r < (chunk[cur_state].size() / conf.row_len); ++r) {
                for(int c = 0; c < conf.row_len; ++c) {
                    int i = r * conf.row_len + c;

                    int neighbours_alive = 0;
                    std::vector<std::array<int, BOARD_DIM>> neighborhood({
                        {-1, -1}, {-1, 0}, {-1, 1},
                        { 0, -1},          { 0, 1},
                        { 1, -1}, { 1, 0}, { 1, 1},
                    });

                    for(const auto &n: neighborhood) {
                        uint8_t s = get_state_wa(chunk[cur_state], neigh_rows, conf.row_len, r + n[0], c + n[1]);
                        //uint8_t s = get_state_closed(rank, size, chunk[cur_state], neigh_rows, conf.row_len, r + n[0], c + n[1]);
                        neighbours_alive += s;
                    }


                    if(chunk[cur_state][i] == 1 && neighbours_alive < 2) {
                        chunk[next_state][i] = 0;
                    }
                    else if(chunk[cur_state][i] == 1 && (neighbours_alive == 2 || neighbours_alive == 3)) {
                        chunk[next_state][i] = 1;
                    }
                    else if(chunk[cur_state][i] == 1 && neighbours_alive > 3) {
                        chunk[next_state][i] = 0;
                    }
                    else if(chunk[cur_state][i] == 0 && neighbours_alive == 3) {
                        chunk[next_state][i] = 1;
                    }

                    if(rank == 1) {
                        LOG(rank, "(%d, %d) %d n=%d prev=%d next=%d", r, c, i, neighbours_alive, chunk[cur_state][i], chunk[next_state][i]);
                    }
                }
            }

            std::swap(next_state, cur_state);
        }
    }

    MPI_Gatherv(chunk[cur_state].data(), proc_cells[rank], MPI_CHAR, board.data(), proc_cells, displs, MPI_CHAR, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        print_board(rank, proc_rows, board, conf, true);
    }

    MPI_Finalize();

    return 0;
}
