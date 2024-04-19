/**
 * @file life.cpp
 * @author Vojtech Dvorak (xdvora3o@fit.vutbr.cz)
 * 
 * ## About
 * 
 * Distributed implementation of Game Of Life in C++ with OpenMPI. Program reads the initial board state from
 * text file and performs N steps (both things are specififed by cmdline arguments). The board is parsed by the
 * root processor, separated to the chunks and distributed (approx. uniformly) over all available
 * processors. Chunks are group of adjacent rows in the board. Communication of processors managing chunks of
 * the board is supplied by cartesian communicator.
 * 
 * By default there are wrap-around edges (but closed edges can be requested by '-c' option). The board is not
 * extensible during the computation - it has dimensions of the initial board from the file. Any rectangular
 * shape of the board is supported (it has not to be a square). 
 * 
 * 
 * ## Usage
 * 
 * ```
 * mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp
 * mpirun --prefix /usr/local/share/OpenMPI -np <NUMBER_OF_PROCS> life <INITIAL_BOARD> <STEPS> [OPTIONS]
 * ```
 * 
 * `NUMBER_OF_PROCS` - Number of processors, all (non-negative) numbers are supported, but situation with 4 procs is the most tested
 * 
 * `INITIAL_BOARD` - Path to the file with initial board
 * 
 * `STEPS` - Number of performed steps
 * 
 * `OPTIONS` (optional) - Flags '-c' and/or '-s':
 *      `-c`  Requests closed edges
 *      `-s` "Silent" mode. Program prints only raw board after the computation (good for testing and restoring the state of the game board from the file) 
 * 
 * 
 * ## Extensions
 * 
 * - Any rectangular board is supported (for example 100x2)
 * - Both edge types are implemented (wrap-around by default, closed can be requested by the flag) 
 * - Works with any number of processor 
 * 
 */

#include <mpi.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>

//#define DEBUG

#define MIN_ARG_NUM 2 //< Amount of mandatory arguments
#define HISTORY_LEN 1 //< Number of steps in the past which will be used for te computation of the next state
#define BOARD_DIM 2 //< Dimensions of game board

typedef enum { BOT, TOP, NEIGH_NUM } neighbour_index_t; //< Indeces for neighbor rows


#ifdef DEBUG
    #define LOG(rank, ...)\
        fprintf(stderr, "Rank: %d  ", rank);\
        fprintf(stderr, __VA_ARGS__);\
        fprintf(stderr, "\n");\
        fflush(stderr);
#else
    #define LOG(...) ""
#endif


/**
 * @brief Coordinates for representing position on the game board
 * 
 */
struct coords_t {
    int row; //< row index
    int col; //< column index

    coords_t operator+(const coords_t &rhs) const {
        return {row + rhs.row, col + rhs.col}; 
    }

};


/**
 * @brief Configuration of an instance of the game of life
 * 
 */
typedef struct game_config {
    size_t row_len; //< Length of each row (only rectangular boards are supported)
    size_t row_num; //< Number of rows
    size_t rows_per_proc; //< Average number of rows per process 
    size_t iteration_n; //< Number of steps of the game
    bool silent; //< Flag for modification of print of the result (can be used for disabling printing of processor indeces)
    uint8_t (*get_state_f)( //< Function for computetion of the next state (wrap-around/closed edges)
        int, int,
        const std::vector<uint8_t> &, const std::vector<uint8_t> *,
        int, const coords_t &);
} game_config_t;


/**
 * @brief Parses initial game board from the text file, only rectangular shape
 * of the game board is supported (row len and row num may be different!)
 * 
 * @param input_file_path Path to the file with the game board
 * @param board Output parameter, the result
 * @param row_len Output param, number of elements in each row of the board
 * @param row_num Output param, number of rows
 * @return int error code, 0 if everything went well
 */
int parse_board(
    const char* input_file_path,
    std::vector<uint8_t> &board,
    size_t &row_len,
    size_t &row_num)
{
    std::ifstream input_file(input_file_path, std::ios::in);
    if(input_file.fail()) {
        std::cerr << "Error: Unable to read file " << input_file_path << "!" << std::endl;
        return 2;
    }

    board.clear();
    row_len = -1;

    char c;
    bool trailing_spaces = false;
    size_t line_i = 0, cur_row_len = 0;
    while(!input_file.eof()) {
        input_file.get(c);
        if(c == '\n' || c == '\r' || input_file.eof()) { //< Line may be ended by LF/CRLF or by EOF (I tolerate missing trailing newline)
            if(c == '\r') {
                input_file.get(c);
            }

            if(line_i == 0) {
                row_len = cur_row_len;
            }
            else {
                if(cur_row_len == 0 && !trailing_spaces) { //< Final line was reached there should be only empty lines from now
                    trailing_spaces = true;
                }

                if(row_len != cur_row_len && !trailing_spaces) {
                    std::cerr << "Error: line " << (line_i + 1);
                    std::cerr << ": The input board must have a rectangular shape!" << std::endl;
                    return 1;
                }
            }

            if(trailing_spaces) {
                continue;
            }

            line_i++;
            cur_row_len = 0;
            if(input_file.eof()) {
                break;
            }
        }
        else {
            board.push_back(c == '1' ? 1 : 0); //< Characters that are not '1' are interpreted as zeros
            cur_row_len++;
            trailing_spaces = false;
        }
    }

    row_num = line_i;
    return 0;
}


/**
 * @brief Prints board to the output
 * 
 * @param size Number of processors
 * @param proc_rows Array with number of rows assigned to the each processor
 * @param board Board to be printed
 * @param config Configuration struct
 */
void print_board(
    int size,
    int *proc_rows,
    std::vector<uint8_t> board,
    game_config_t &config)
{
    size_t p = 0;
    for(size_t i = 0, proc_row_cnt = 0; i < config.row_num; ++i, ++proc_row_cnt) {
        if(!config.silent) { //< Print prefixes with processor ranks only if silent mode is deactivated
            while(proc_rows[p] <= proc_row_cnt) {
                proc_row_cnt = 0;
                if(p + 1 > size) {
                    break;
                }

                ++p;
            }

            std::cout << p << ": ";
        }

        // Print every column in the current row
        for(size_t j = 0; j < config.row_len; ++j) {
            std::cout << (int)board[i * config.row_len + j];
        }

        std::cout << std::endl;
    }
}


/**
 * @brief Gets state of one cell, while considering wrap-around edges
 * NOTE: It presumes, that caller will ask for state of cell corresponding
 * processor OR the celll will be available in neigh_rows params
 * 
 * @param rank Rank of current processor
 * @param size The number of processors
 * @param chunk Chunk of the board dedicated to the current processor
 * @param neigh_rows Vectors with adjacent rows
 * @param row_len Length of each row
 * @param orig_coords Coordinates of the cell in !the chunk!
 * @return uint8_t State of the cell
 */
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
    if(coords.col >= row_len) { //< If col is greater or eq than number of cells in rows make wrap-around 
        coords.col = 0;
    }
    else if(coords.col < 0) {
        coords.col += row_len; //< Same if there is negative column coordinate
    }

    // With rows is situation a little bit complicated - we must borrow state of the cells from the neighbors
    const uint8_t *source = chunk.data();
    if(coords.row >= row_num) { // If row coord is greater or eq than number of rows in chunk look at the adjacent row
        coords.row = 0;
        source = &(neigh_rows[TOP][row_len]);
    }
    else if(coords.row < 0) { // Similarly when coord is negative
        coords.row = 0;
        source = &(neigh_rows[BOT][0]);
    }

    return source[coords.row * row_len + coords.col];
}


/**
 * @brief Gets state of one cell, while considering closed edges (there
 * are 0 outside the board)
 * NOTE: It presumes, that caller will ask for state of cell corresponding
 * processor OR the celll will be available in neigh_rows params
 * 
 * @param rank Rank of current processor
 * @param size The number of processors
 * @param chunk Chunk of the board managed by the current processor
 * @param neigh_rows Vectors with adjacent rows
 * @param row_len Length of each row (number of elements in each row)
 * @param orig_coords Coordinates of the cell in !the chunk!
 * @return uint8_t The value of the cell
 */
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
    if(coords.col >= row_len || coords.col < 0) { // Cell is outside of the board -> 0
        return 0;
    }

    const uint8_t *source = chunk.data();
    if(coords.row >= row_num) { // Cell is outside of the chunk
        if(rank >= size - 1) { // Cell is outside of the board -> 0
            return 0;
        }

        coords.row = 0;
        source = &(neigh_rows[TOP][row_len]);
    }
    else if(coords.row < 0) { // The cell is outside of the chunk
        if(rank == 0) { // Cell is outside of the board -> 0
            return 0;
        }

        coords.row = 0;
        source = &(neigh_rows[BOT][0]);
    }

    return source[coords.row * row_len + coords.col];
}


/**
 * @brief Computes the next state of the cell specified by the coords
 * 
 * @param rank Rank of the current processor
 * @param size Total number of the processors
 * @param cur_state_chunk Chunk managed by the current processor with current state
 * @param neigh_rows Adjacent rows of the chunk
 * @param config Game configuration structure
 * @param coords Coordinates of the cell its nest state we want to compute
 * @return uint8_t The next value of the cell
 */
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

    // Relative coordinates of cells in the neighborhood
    std::vector<coords_t> neighborhood({
        {-1, -1}, {-1, 0}, {-1, 1},
        { 0, -1},          { 0, 1},
        { 1, -1}, { 1, 0}, { 1, 1},
    });

    // Collect states of all cells in the neighborhood
    int neighbors_alive = 0;
    for(const auto &n: neighborhood) {
        coords_t neigbor_coords = coords + n;
        uint8_t s = config.get_state_f(
            rank, size, cur_state_chunk, neigh_rows, config.row_len, neigbor_coords
        );
        neighbors_alive += s;
    }

    // Rules of the game of live (specified in assignment)
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


/**
 * @brief Compute the state of the cells in the chunk after number of steps
 * thatis specified in the config structure
 * 
 * @param rank Rank of the current processor
 * @param size Total number of processors
 * @param chunk Chunks managed by the current processors (current state chunk + history)
 * @param board_comm Cartesian communicator with all cells in the board
 * @param config Configuration structure
 * @return int Index of the chunk with the current state at the end of the game
 */
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

    // Every procesor uses array of chunks to easily acces the current state of cells and stores their next state
    int cur_state = 0, next_state = 1; //< There are only two chunks (one with current state and one with the next state)

    // Prepare vectors for storing adjacent rows of the chunk - we will just need their state, their next state will be computed by another processor
    std::vector<uint8_t> board, neigh_rows[NEIGH_NUM];
    neigh_rows[BOT].resize(NEIGH_NUM * config.row_len);
    neigh_rows[TOP].resize(NEIGH_NUM * config.row_len);

    for(size_t i = 0; i < config.iteration_n; ++i) {
        uint8_t* last_row_ptr = &(chunk[cur_state].data()[chunk[cur_state].size() - config.row_len]);

        // Fetch adjacent rows (! we will get edges rows from the both neighbors)
        MPI_Neighbor_allgather(last_row_ptr, config.row_len, MPI_CHAR, neigh_rows[BOT].data(), config.row_len, MPI_CHAR, board_comm);
        MPI_Neighbor_allgather(chunk[cur_state].data(), config.row_len, MPI_CHAR, neigh_rows[TOP].data(), config.row_len, MPI_CHAR, board_comm);

        memset(chunk[next_state].data(), 0, chunk[next_state].size()); //< Reset chunk with next state 

        // Iterate over all cells in the chunk
        for(int r = 0; r < (chunk[cur_state].size() / config.row_len); ++r) {
            for(int c = 0; c < config.row_len; ++c) {
                coords_t coords = {r, c};

                // Compute next state of the current cell
                chunk[next_state][r * config.row_len + c] = get_next_state(rank, size, chunk[cur_state], neigh_rows, config, coords);
            }
        }

        std::swap(next_state, cur_state);
    }

    return cur_state;
}


/**
 * @brief Creates board ccooomunicator - cartesian communicator with only those
 * processors that manage any rows
 * 
 * @param rank Rank of the current processor
 * @param proc_rows Array with numbers of rows for each processor
 * @return MPI_Comm Board communicator
 */
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


/**
 * @brief Gets details about chunk of the current processor
 * 
 * @param size Total number of processors
 * @param displs Displacement of chunk cells in the whole board (output)
 * @param proc_cells Arra with umber of cells managed by each processor
 * (output)
 * @param proc_rows Number of rows managed by the current processor
 * @param config Configuration structure, that stores detail about 
 * the current game instace
 */
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


/**
 * @brief Scatter chunks - send chunks of the board (some amount of rows) 
 * to each processor, in the world
 * 
 * @param rank Rank of the currenkt processor
 * @param board Board to be distributed
 * @param chunk Chunks of the board managed by the current processor (output)
 * @param displs Displacement of chunk cells in the whole board (output)
 * @param proc_cells Number of cells managed by each processor (output)
 */
void scatter_chunks(
    int rank,
    const std::vector<uint8_t> &board,
    std::vector<uint8_t> *chunk,
    const int *displs,
    const int *proc_cells)
{
    //LOG(rank, "%d\n", proc_cells[rank]);
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


/**
 * @brief Gets the all chunks managed by all processors 
 * 
 * @param rank Rank of the curren proccesors
 * @param board Board which will be filled by the chunks (output)
 * @param end_state_chunk Chunk with the final state managed by current processor 
 * @param displs Displacement of each chunk in the whole board
 * @param proc_cells Array with number of cells managed by each processors 
 */
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


/**
 * @brief Prepares game of life (should be called only by the one processor),
 * parses board in the file and computes other configuration of the board
 * 
 * @param argc Argc
 * @param argv Argv
 * @param size The total number of processors
 * @param board Game board filled by the state from the file (output)
 * @param proc_rows Array with number of rows managed by each processor (output)
 * @param config Configuration structure of the current game instance
 * @return int 0 if everything went well, non-zero if there were any error
 */
int prepare_game(
    int argc,
    char **argv,
    int size,
    std::vector<uint8_t> &board,
    int *proc_rows,
    game_config_t &config)
{
    // Parse cmdline args
    if(argc < MIN_ARG_NUM + 1) {
        std::cerr << "Usage error!" << std::endl;
        std::cerr << "USAGE: ./life <path-to-board> <it-num> [OPTIONS]" << std::endl;
        return 1;
    }

    const char *board_path = argv[1];

    try {
        config.iteration_n = std::stoi(argv[2]);
    }
    catch(const std::exception &e) {
        std::cerr << "Error: Bad iteration number " << argv[2] << std::endl;
        return 1;
    }

    // Parse initial game board
    int ret = parse_board(argv[1], board, config.row_len, config.row_num);
    if(ret) {
        return ret;
    }

    // Compute number of rows for each processor
    if(size == 1) {
        config.rows_per_proc = config.row_num;
        proc_rows[0] = config.rows_per_proc;
    }
    else if(size == 2) {
        config.rows_per_proc = config.row_num / 2;
        proc_rows[0] = config.rows_per_proc;
        proc_rows[1] = config.row_num - config.rows_per_proc;
    }
    else {
        size_t remaining_rows = config.row_num; 
        config.rows_per_proc = config.row_num / (size - 1);
        for(size_t i = 0; remaining_rows >= config.rows_per_proc && i < size - 1; ++i) {
            proc_rows[i] = config.rows_per_proc;
            remaining_rows -= config.rows_per_proc;
        }
        proc_rows[size - 1] = remaining_rows;
    }

    LOG(0, "Board size: %ldx%ld", config.row_len, config.row_num);
    LOG(0, "Rows per process: %ld", config.rows_per_proc);
    LOG(0, "Remaining rows: %d", proc_rows[size - 1]);

    return 0;
}


/**
 * @brief Configurates each processors
 * 
 * @param argc Argc 
 * @param argv Argv
 * @param size The total number of processors
 * @param proc_rows Number of rows managed by each processor
 * @param config Configuration of current game of life instace 
 */
void configurate_procs(
    int argc,
    char **argv,
    int size,
    int *proc_rows,
    game_config_t *config)
{
    MPI_Bcast(proc_rows, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(config, sizeof(game_config_t), MPI_BYTE, 0, MPI_COMM_WORLD);

    config->get_state_f = get_state_wa;
    config->silent = false;
    for(int i = MIN_ARG_NUM + 1; i < argc; ++i) {
        if(!strcmp(argv[i], "-c")) { //< There is an extension - user can choose between wrap-around (default) and closed edged by cmdline argument -c
            config->get_state_f = get_state_closed;
        }
        else if(!strcmp(argv[i], "-s")) {  //< There is an extension - user can deactivate processor labels from the output of the program by -s c
            config->silent = true;
        }
    }
}



int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size;

    // Get rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<uint8_t> board, chunk[HISTORY_LEN + 1]; // Complete board and chunk (chunks are parts of the board each managed by different processors)
    int proc_rows[size], result = 0;
    game_config_t config;
    if(rank == 0) { // The root processor
        result = prepare_game(argc, argv, size, board, proc_rows, config);
    }

    MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the result of the initilaization, end program if there were any error
    if(result) {
        MPI_Finalize();
        return result;
    }

    configurate_procs(argc, argv, size, proc_rows, &config); // Configure each processor

    int displs[size], proc_cells[size];
    get_chunk_specs(size, displs, proc_cells, proc_rows, config); // Get details of chunks (compute number of cells from number of rows etc.)

    scatter_chunks(rank, board, chunk, displs, proc_cells); // Distribute chunks of the board (rank 0 disitributes chunks to other processors)

    MPI_Comm board_comm = create_board_comm(rank, proc_rows); // Create cartesian communicator for communcation of processors with some chunk

    int cur_state = compute(rank, size, chunk, board_comm, config);

    gather_chunks(rank, board, chunk[cur_state], displs, proc_cells); // Get all chunks back after the game and recostruct the board

    if(rank == 0) { // Rank 0 have now complete board and it is able to print it
        print_board(size, proc_rows, board, config);
    }

    MPI_Finalize();

    return 0;
}

/***                          End of the life.cpp                              */
