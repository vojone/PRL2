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

//#define DEBUG

#define MIN_ARG_NUM 2
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


typedef struct game_confing {
    size_t row_len;
    size_t row_num;
    size_t rows_per_proc;
    size_t iteration_n;
} game_config_t;


struct coords_t {
    int row;
    int col;

    coords_t operator+(const coords_t &rhs) const {
        return {row + rhs.row, col + rhs.col}; 
    }

};


int parse_board(
    const char* input_file_path,
    std::vector<uint8_t> &board,
    size_t &row_len,
    size_t &row_num)
{
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


void print_board(
    int rank,
    int *proc_rows,
    std::vector<uint8_t> board,
    game_config_t &config,
    bool print_prefixes)
{
    size_t p = 0;
    for(size_t i = 0, proc_row_cnt = 0; i < config.row_num; ++i, ++proc_row_cnt) {
        if(print_prefixes) {
            while(proc_rows[p] <= proc_row_cnt) {
                proc_row_cnt = 0;
                ++p;
            }

            std::cout << p << ": ";
        }

        for(size_t j = 0; j < config.row_len; ++j) {
            std::cout << (int)board[i * config.row_len + j];
        }

        std::cout << std::endl;
    }
}


uint8_t get_state_wa(
    int rank,
    int size,
    const std::vector<uint8_t> &chunk,
    const std::vector<uint8_t> *neigh_rows,
    int row_len,
    const coords_t &orig_coords)
{
    coords_t coords = orig_coords;

    int row_num = chunk.size() / row_len;
    if(coords.col >= row_len) {
        coords.col = 0;
    }
    else if(coords.col < 0) {
        coords.col += row_len;
    }

    const uint8_t *source = chunk.data();
    if(coords.row >= row_num) {
        coords.row = 0;
        source = &(neigh_rows[TOP][row_len]);
    }
    else if(coords.row < 0) {
        coords.row = 0;
        source = &(neigh_rows[BOT][0]);
    }

    return source[coords.row * row_len + coords.col];
}


uint8_t get_state_closed(
    int rank,
    int size,
    const std::vector<uint8_t> &chunk,
    const std::vector<uint8_t> *neigh_rows,
    int row_len,
    const coords_t &orig_coords)
{
    coords_t coords = orig_coords;
    int row_num = chunk.size() / row_len;
    if(coords.col >= row_len || coords.col < 0) {
        return 0;
    }

    const uint8_t *source = chunk.data();
    if(coords.row >= row_num) {
        if(rank >= size - 1) {
            return 0;
        }

        coords.row = 0;
        source = &(neigh_rows[TOP][row_len]);
    }
    else if(coords.row < 0) {
        if(rank == 0) {
            return 0;
        }

        coords.row = 0;
        source = &(neigh_rows[BOT][0]);
    }

    return source[coords.row * row_len + coords.col];
}


uint8_t get_next_state(
    int rank,
    int size,
    const std::vector<uint8_t> &cur_state_chunk,
    const std::vector<uint8_t> *neigh_rows,
    const game_config_t &config,
    const coords_t &coords) 
{
    int i = coords.row * config.row_len + coords.col;
    uint8_t result = cur_state_chunk[i];

    std::vector<coords_t> neighborhood({
        {-1, -1}, {-1, 0}, {-1, 1},
        { 0, -1},          { 0, 1},
        { 1, -1}, { 1, 0}, { 1, 1},
    });

    int neighbors_alive = 0;
    for(const auto &n: neighborhood) {
        coords_t neigbor_coords = coords + n;
        //uint8_t s = get_state_wa(rank, size, cur_state_chunk, neigh_rows, config.row_len, neigbor_coords);
        uint8_t s = get_state_closed(rank, size, cur_state_chunk, neigh_rows, config.row_len, neigbor_coords);
        neighbors_alive += s;
    }

    if(cur_state_chunk[i] == 1 && neighbors_alive < 2) {
        result = 0;
    }
    else if(cur_state_chunk[i] == 1 && (neighbors_alive == 2 || neighbors_alive == 3)) {
        result = 1;
    }
    else if(cur_state_chunk[i] == 1 && neighbors_alive > 3) {
        result = 0;
    }
    else if(cur_state_chunk[i] == 0 && neighbors_alive == 3) {
        result = 1;
    }

    return result;
}


int compute(
    int rank,
    int size,
    std::vector<uint8_t> *chunk,
    MPI_Comm board_comm,
    const game_config_t &config) 
{
    if(chunk[0].size() == 0) {
        return 0;
    }

    int cur_state = 0, next_state = 1;

    std::vector<uint8_t> board, neigh_rows[NEIGH_NUM];
    neigh_rows[BOT].resize(NEIGH_NUM * config.row_len);
    neigh_rows[TOP].resize(NEIGH_NUM * config.row_len);

    for(size_t i = 0; i < config.iteration_n; ++i) {
        uint8_t* last_row_ptr = &(chunk[cur_state].data()[chunk[cur_state].size() - config.row_len]);

        MPI_Neighbor_allgather(last_row_ptr, config.row_len, MPI_CHAR, neigh_rows[BOT].data(), config.row_len, MPI_CHAR, board_comm);
        MPI_Neighbor_allgather(chunk[cur_state].data(), config.row_len, MPI_CHAR, neigh_rows[TOP].data(), config.row_len, MPI_CHAR, board_comm);

        memset(chunk[next_state].data(), 0, chunk[next_state].size());
        for(int r = 0; r < (chunk[cur_state].size() / config.row_len); ++r) {
            for(int c = 0; c < config.row_len; ++c) {
                coords_t coords = {r, c};
                chunk[next_state][r * config.row_len + c] = get_next_state(rank, size, chunk[cur_state], neigh_rows, config, coords);
            }
        }

        std::swap(next_state, cur_state);
    }

    return cur_state;
}


MPI_Comm create_board_comm(int rank, const int *proc_rows)
{
    MPI_Comm chunk_holder_comm, board_comm;
    MPI_Comm_split(MPI_COMM_WORLD, proc_rows[rank] > 0, rank, &chunk_holder_comm);

    int chunk_holder_num;
    MPI_Comm_size(chunk_holder_comm, &chunk_holder_num);

    int dims[BOARD_DIM] = {chunk_holder_num}, periodic[1] = {1};
    MPI_Cart_create(chunk_holder_comm, 1, dims, periodic, 0, &board_comm);


    return board_comm;
}


void get_chunk_specs(
    int size,
    int *displs,
    int *proc_cells,
    const int *proc_rows,
    const game_config_t &config)
{
    for(size_t i = 0; i < size; ++i) {
        proc_cells[i] = config.row_len * proc_rows[i];
        displs[i] = i * config.row_len * config.rows_per_proc;
    }
}


void scatter_chunks(
    int rank,
    const std::vector<uint8_t> &board,
    std::vector<uint8_t> *chunk,
    const int *displs,
    const int *proc_cells)
{
    chunk[0].resize(proc_cells[rank]);
    
    MPI_Scatterv(
        board.data(),
        proc_cells,
        displs,
        MPI_CHAR,
        chunk[0].data(),
        proc_cells[rank],
        MPI_CHAR,
        0,
        MPI_COMM_WORLD
    );
    chunk[1].assign(chunk[0].begin(), chunk[0].end());
}


void gather_chunks(
    int rank,
    std::vector<uint8_t> &board,
    const std::vector<uint8_t> &end_state_chunk,
    const int *displs,
    const int *proc_cells)
{
    MPI_Gatherv(
        end_state_chunk.data(),
        proc_cells[rank],MPI_CHAR,
        board.data(),
        proc_cells,
        displs,
        MPI_CHAR,
        0,
        MPI_COMM_WORLD
    );
}


int prepare_game(
    int argc,
    char **argv,
    int size,
    std::vector<uint8_t> &board,
    int *proc_rows,
    game_config_t &config)
{
    if(argc < MIN_ARG_NUM + 1) {
        std::cerr << "USAGE: ./life <path-to-board> <it-num>" << std::endl;
        return 2;
    }

    const char *board_path = argv[1];
    config.iteration_n = std::stoi(argv[2]);


    int ret = parse_board(argv[1], board, config.row_len, config.row_num);
    if(ret) {
        MPI_Finalize();
        return ret;
    }

    config.rows_per_proc = config.row_num / (size - 1);
    size_t remaining_rows = config.row_num; 
    for(size_t i = 0; remaining_rows >= config.rows_per_proc && i < size - 1; ++i) {
        proc_rows[i] = config.rows_per_proc;
        remaining_rows -= config.rows_per_proc;
    }

    proc_rows[size - 1] = remaining_rows;

    LOG(0, "Board size: %ldx%ld", config.row_len, config.row_num);
    LOG(0, "Rows per process: %ld", config.rows_per_proc);
    LOG(0, "Remaining rows: %d", proc_rows[size - 1]);

    return 0;
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size;

    // Get rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<uint8_t> board, chunk[HISTORY_LEN];
    int proc_rows[size];
    game_config_t config;
    if(rank == 0) { // The root processor processor
        int result = prepare_game(argc, argv, size, board, proc_rows, config);
        if(result) {
            MPI_Finalize();
            return result;
        }
    }

    MPI_Bcast(proc_rows, size, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&config, sizeof(config), MPI_BYTE, 0, MPI_COMM_WORLD);

    int displs[size], proc_cells[size];
    get_chunk_specs(size, displs, proc_cells, proc_rows, config);

    scatter_chunks(rank, board, chunk, displs, proc_cells);


    MPI_Comm board_comm = create_board_comm(rank, proc_rows);

    int cur_state = compute(rank, size, chunk, board_comm, config);

    gather_chunks(rank, board, chunk[cur_state], displs, proc_cells);

    if(rank == 0) {
        print_board(rank, proc_rows, board, config, true);
    }

    MPI_Finalize();

    return 0;
}
