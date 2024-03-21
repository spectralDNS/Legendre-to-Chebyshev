#!/bin/bash

# Loop through input arguments n=2 to n=10
output=$(l2c -N256 -d0 -r10000 -s32 -l1 -v0)
number=$(echo "$output" | awk '{print $10}')
echo "[256, $number],"
sleep 1

for (( n=1; n<8; n++ ))
do
    result=$(echo "64*2^($n+2)" | bc)
    output=$(l2c -N$result -d0 -r10000 -s64 -l1 -v0)
    number=$(echo "$output" | awk '{print $10}')
    echo "[$result, $number],"
    sleep 1
done

for (( n=8; n<12; n++ ))
do
    result=$(echo "64*2^($n+2)" | bc)
    output=$(l2c -N$result -d0 -r1000 -s64 -l1 -v0)
    number=$(echo "$output" | awk '{print $10}')
    echo "[$result, $number],"
    sleep 1
done

for (( n=12; n<=15; n++ ))
do
    result=$(echo "64*2^($n+2)" | bc)
    output=$(l2c -N$result -d0 -r100 -s64 -l1 -v0)
    number=$(echo "$output" | awk '{print $10}')
    echo "[$result, $number],"
    sleep 1
done