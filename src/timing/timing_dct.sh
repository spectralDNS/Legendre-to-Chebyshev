#!/bin/bash

for (( n=0; n<=15; n++ ))
do
    result=$(echo "64*2^($n+0)" | bc)
    output=$(l2c -N$result -d7 -r100)
    number=$(echo "$output" | awk '{print $2}')
    echo "[$result, $number],"
    sleep 1
done
