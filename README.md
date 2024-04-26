# PRL - 2. project - Game of Life

## Author

Vojtěch Dvořák (xdvora3o)

## About 
Distributed implementation of Game Of Life in C++ with OpenMPI. Program reads the initial board state from
text file and performs N steps (both things are specififed by cmdline arguments). The board is parsed by the
root processor, separated to the chunks and distributed (approx. uniformly) over all available
processors. Chunks are group of adjacent rows in the board. Communication of processors managing chunks of
the board is supplied by cartesian communicator.
 
By default there are wrap-around edges (but closed edges can be requested by '-c' option). The board is not
extensible during the computation - it has dimensions of the initial board from the file. Any rectangular
shape of the board is supported (it has not to be a square). 


## Usage
 
 ```
 mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp
 mpirun --prefix /usr/local/share/OpenMPI -np <NUMBER_OF_PROCS> life <INITIAL_BOARD> <STEPS> [OPTIONS]
 ```
 
 `NUMBER_OF_PROCS` - Number of processors, all (non-negative) numbers are supported, but situation with 4 procs is the most tested
 
 `INITIAL_BOARD` - Path to the file with initial board
 
 `STEPS` - Number of performed steps
 
 `OPTIONS` (optional) - Flags '-c' and/or '-s':
      
      `-c`  Requests closed edges
      
      `-s` "Silent" mode. Program prints only raw board after the computation (good for testing and restoring the state of the game board from the file) 
      
 
## Extensions
 
- Any rectangular board is supported (for example 100x2)

- Both edge types are implemented (wrap-around by default, closed can be requested by the flag) 

- Works with any number of processor 
  
