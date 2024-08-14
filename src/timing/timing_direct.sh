#!/bin/bash

for (( n=1; n<10; n++ ))
do
    result=$(echo "32*2^($n)" | bc)
    output=$(l2c -N$result -d6 -r10000 -v0)
    number=$(echo "$output" | awk '{print $10}')
    echo "[$result, $number],"
    sleep 1
done