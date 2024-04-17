#!/bin/bash

if [ $# -lt 2 ]; then
    echo "USAGE: ./test.sh <initial-board-path> <iteration-num> <extra-args>"
    exit 1
fi

#preklad zdrojoveho souboru
mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp

#spusteni programu
# shellcheck disable=SC2086
mpirun --prefix /usr/local/share/OpenMPI  -np 4 life "$1" "$2" $3

#uklid
rm -f life
