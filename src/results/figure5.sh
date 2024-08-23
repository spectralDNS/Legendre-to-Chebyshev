#!/bin/bash

# Run this file to produce Fig. 5 in the paper 
# Set any of the reruns below to yes in order to rerun and recreate that dataset
# Otherwise, use already computed data stored below in ascii.
# This bash-file generates a python-file with speed-data and plots fig 5
# and stores it as threadscaling.png. If threadscaling.png exists it will
# be overwritten.
# As long as the trusted data are stored in this file it is perfectly safe
# to rerun either the fmm or the dct. The new data will only appear in a 
# new generated python file (see filename below).

rerun_fmm="no"
rerun_dct="no"
smax=32
model=0

filename="speedthreads_USEACC_${USE_ACCELERATE}_${RANDOM}_s${smax}_model${model}_rerunfmm_${rerun_fmm}_rerundct_${rerun_dct}.py"
echo $filename
echo "import numpy as np" > $filename
echo import matplotlib.pyplot as plt >> $filename
if [ "$rerun_fmm" = "yes" ]; then
for (( threads=1; threads<16; threads=threads*2 ))
do
echo "t${threads} = np.array([" >> $filename
for (( n=8; n<=23; n++ ))
do
    result=$(echo "2^($n)" | bc)
    repeat=$(echo "2^(26-$n)" | bc)
    output=$(l2c -N$result -d0 -r$repeat -s$smax -l$model -v0 -t$threads)
    number=$(echo "$output" | awk '{print $10}')
    echo "[$result, $number]," >> $filename
    sleep 0.1
done
echo "])" >> $filename
done
else 
cat << EOF >> $filename 
t1 = np.array([
  [256, 2.5432e-05],
[512, 5.6912e-05],
[1024, 1.0378e-04],
[2048, 1.9178e-04],
[4096, 3.5994e-04],
[8192, 7.2041e-04],
[16384, 1.5880e-03],
[32768, 2.7371e-03],
[65536, 5.2833e-03],
[131072, 1.0286e-02],
[262144, 2.2188e-02],
[524288, 4.5583e-02],
[1048576, 9.1787e-02],
[2097152, 1.8770e-01],
[4194304, 3.7999e-01],
[8388608, 7.6695e-01]
])
t2 = np.array([
  [256, 3.5100e-05],
[512, 6.7710e-05],
[1024, 1.0886e-04],
[2048, 1.7253e-04],
[4096, 2.7519e-04],
[8192, 4.7911e-04],
[16384, 9.2138e-04],
[32768, 1.7081e-03],
[65536, 3.2525e-03],
[131072, 6.2187e-03],
[262144, 1.2646e-02],
[524288, 2.5416e-02],
[1048576, 5.1586e-02],
[2097152, 1.0628e-01],
[4194304, 2.1335e-01],
[8388608, 4.2672e-01]
])
t4 = np.array([
[256, 3.4232e-05],
[512, 6.3518e-05],
[1024, 9.8978e-05],
[2048, 1.4294e-04],
[4096, 2.1336e-04],
[8192, 3.3536e-04],
[16384, 5.6715e-04],
[32768, 1.0185e-03],
[65536, 1.8988e-03],
[131072, 3.6117e-03],
[262144, 7.2761e-03],
[524288, 1.4235e-02],
[1048576, 2.9020e-02],
[2097152, 5.8746e-02],
[4194304, 1.1812e-01],
[8388608, 2.3558e-01]
])
t8 = np.array([
[256, 4.0460e-05],
[512, 6.9500e-05],
[1024, 1.0072e-04],
[2048, 1.3552e-04],
[4096, 1.9006e-04],
[8192, 2.7106e-04],
[16384, 4.2366e-04],
[32768, 7.0978e-04],
[65536, 1.2622e-03],
[131072, 2.2207e-03],
[262144, 4.4938e-03],
[524288, 8.4946e-03],
[1048576, 1.6889e-02],
[2097152, 3.4687e-02],
[4194304, 6.9643e-02],
[8388608, 1.4133e-01]
])
EOF
fi

# dct
if [ "$rerun_dct" = "yes" ]; then
for (( threads=1; threads<16; threads=threads*2 ))
do
echo "dct${threads} = np.array([" >> $filename
for (( n=8; n<=23; n++ ))
do
    result=$(echo "2^($n)" | bc)
    repeat=$(echo "2^(26-$n)" | bc)
    output=$(l2c -N$result -d7 -r$repeat -t$threads)
    number=$(echo "$output" | awk '{print $11}')
    echo "[$result, $number]," >> $filename
    sleep 0.1
done
echo "])" >> $filename
done
else 
cat << EOF >> $filename
dct1 = np.array([
[256, 2.212000e-06],
[512, 4.344000e-06],
[1024, 7.060000e-06],
[2048, 1.506400e-05],
[4096, 3.368800e-05],
[8192, 7.384000e-05],
[16384, 1.590400e-04],
[32768, 5.832000e-04],
[65536, 7.917600e-04],
[131072, 1.753400e-03],
[262144, 3.801908e-03],
[524288, 8.684952e-03],
[1048576, 2.009936e-02],
[2097152, 4.412776e-02],
[4194304, 1.216114e-01],
[8388608, 2.630520e-01]
])
dct2 = np.array([
[256, 1.546800e-05],
[512, 1.761200e-05],
[1024, 2.140400e-05],
[2048, 2.822000e-05],
[4096, 4.237200e-05],
[8192, 9.172400e-05],
[16384, 2.037280e-04],
[32768, 4.296720e-04],
[65536, 5.720800e-04],
[131072, 1.072036e-03],
[262144, 2.129456e-03],
[524288, 4.875396e-03],
[1048576, 1.075306e-02],
[2097152, 2.390993e-02],
[4194304, 5.948916e-02],
[8388608, 1.281804e-01]
])
dct4 = np.array([
[256, 1.328400e-05],
[512, 1.910400e-05],
[1024, 2.355600e-05],
[2048, 2.970400e-05],
[4096, 4.274800e-05],
[8192, 6.980000e-05],
[16384, 1.519640e-04],
[32768, 3.421120e-04],
[65536, 4.288920e-04],
[131072, 8.241960e-04],
[262144, 1.622348e-03],
[524288, 3.474832e-03],
[1048576, 7.602892e-03],
[2097152, 1.608698e-02],
[4194304, 4.191645e-02],
[8388608, 8.676902e-02]
])
dct8 = np.array([
[256, 2.496800e-05],
[512, 2.922800e-05],
[1024, 2.929600e-05],
[2048, 4.000800e-05],
[4096, 4.983600e-05],
[8192, 7.923200e-05],
[16384, 1.689400e-04],
[32768, 3.097960e-04],
[65536, 4.064040e-04],
[131072, 7.092240e-04],
[262144, 1.405048e-03],
[524288, 3.053396e-03],
[1048576, 6.588604e-03],
[2097152, 1.334960e-02],
[4194304, 3.384324e-02],
[8388608, 6.670523e-02]
])
EOF
fi

cat << EOF >> $filename
fig = plt.figure(figsize=(6, 4))
ax1 = fig.gca()
ax1.semilogx(t1[:, 0], t1[:, 1]/t2[:, 1], 'k')
ax1.semilogx(t1[:, 0], t1[:, 1]/t4[:, 1], 'k:')
ax1.semilogx(t1[:, 0], t1[:, 1]/t8[:, 1], 'k--')
ax1.semilogx(t1[:, 0], dct1[:, 1]/dct2[:, 1], 'ks', ms=4)
ax1.semilogx(t1[:, 0], dct1[:, 1]/dct4[:, 1], 'ko', ms=4)
ax1.semilogx(t1[:, 0], dct1[:, 1]/dct8[:, 1], 'k^', ms=4)
ax1.set_xticks(t1[1::2, 0], labels=np.log2(t1[1::2, 0]).astype(int))
ax1.legend(['2 cores (L2C)', '4 cores (L2C)', '8 cores (L2C)', 
            '2 cores (DCT)', '4 cores (DCT)', '8 cores (DCT)'])
ax1.set_ylabel('Speedup')
ax1.set_xlabel('N')
plt.savefig('threadscaling.png')
plt.show()
EOF

python $filename