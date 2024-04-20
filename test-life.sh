#!/bin/bash

# Testing script for life.cpp (needs OpenMPI to be installed)

shopt -s extglob

TEST_DIR="$1"
EXTRA_ARGS=""
NUMBER_OF_PROC=4

if [ $# -lt 1 ]; then
    printf "USAGE: ./test-life.sh <test-dir> [<number-of-procs>] [<extra-args>]\n"
    printf "\n"
    printf "test-dir\tPath to the directory with test cases\n"
    printf "number-of-procs\tNumber of procs (4 by default)\n"
    printf "extra-args\tExtra arguments for the life binary\n"
    printf "\n"
    exit 1
fi

if [ $# -gt 1 ] && [ ! -z "$2" ]; then
    NUMBER_OF_PROC=$2
fi

if [ $# -gt 2 ]; then
    EXTRA_ARGS=$3
fi

mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp


echo "Extra args: $EXTRA_ARGS"
echo "Proc num: $NUMBER_OF_PROC"

CNT=0
SUCCESS_CNT=0
for CASE in $TEST_DIR/*; do
    echo -n "$CASE: "

    INPUT="$CASE/$(ls "$CASE" | grep ".*.in")"
    TARGET="$CASE/$(ls "$CASE" | grep ".*.out")"
    ITERATION_NUM_FILE="$CASE/$(ls "$CASE" | grep ".*.it")"

    if [ ! -f "$INPUT" ] || [ ! -f "$TARGET" ] || [ ! -f "$ITERATION_NUM_FILE" ]; then
        echo "Incomplete test case!"
    else
        ITERATION_NUM=$(tr -d '[:space:]' < "$ITERATION_NUM_FILE")
        # echo "$ITERATION_NUM_FILE IT: $ITERATION_NUM"
        # shellcheck disable=SC2086
        RESULT=$(mpirun --prefix /usr/local/share/OpenMPI  -np "$NUMBER_OF_PROC" life "$INPUT" "$ITERATION_NUM" $EXTRA_ARGS)
        DIFFS=$(diff -Z "$TARGET" <(echo "$RESULT"))

        if [ -z "$DIFFS" ]; then
            printf "OK"
            SUCCESS_CNT=$((SUCCESS_CNT + 1))
        else
            printf "FAIL\nExpected:\n%s\nGot:\n%s\n" "$(cat "$TARGET")" "$RESULT"
        fi
    fi

    CNT=$((CNT + 1))

    echo
done

echo "Result: $SUCCESS_CNT/$CNT"

rm -f life
