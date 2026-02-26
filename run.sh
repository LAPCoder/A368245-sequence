#!/bin/bash
set -e
# This script is used to verify the numbers found but is it quite slow
# You should use it only if you need, but the C++ code is reliable

# Create buffer in RAM
echo "Creating buffer in /dev/shm"
BUFFER="/dev/shm/data_buffer_$$.bin"
trap "rm -f $BUFFER" EXIT
echo "Using RAM buffer: $BUFFER"

echo "Running C++ and Python in parallel"
./sequence_SIMD 5000 300 > out_hex.pgm 2> "$BUFFER" &
CPPID=$!

while ! test -f "$BUFFER"; do
	sleep 0.001
	echo "Still waiting for CPP to start writting..."
done

echo "Running Python"
python -u checker.py --base 16 < "$BUFFER" > number_found.log &
PYID=$!

wait $CPPID
echo "C++ finished"
wait $PYID
echo "Python finished"
echo "Done!"