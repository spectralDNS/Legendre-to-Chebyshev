#!/bin/bash

model=1
s=32

output=$(l2c -N256 -d9 -r10000 -s32 -l$model -v0)
number=$(echo "$output" | awk '{print $11}')
echo "[256, $number],"
sleep 1

for (( n=1; n<8; n++ ))
do
    result=$(echo "64*2^($n+2)" | bc)
    output=$(l2c -N$result -d9 -r1000 -s$s -l$model -v0)
    number=$(echo "$output" | awk '{print $11}')
    echo "[$result, $number],"
    sleep 1
done

for (( n=8; n<12; n++ ))
do
    result=$(echo "64*2^($n+2)" | bc)
    output=$(l2c -N$result -d9 -r100 -s$s -l$model -v0)
    number=$(echo "$output" | awk '{print $11}')
    echo "[$result, $number],"
    sleep 1
done

for (( n=12; n<=15; n++ ))
do
    result=$(echo "64*2^($n+2)" | bc)
    output=$(l2c -N$result -d9 -r10 -s$s -l$model -v0)
    number=$(echo "$output" | awk '{print $11}')
    echo "[$result, $number],"
    sleep 1
done