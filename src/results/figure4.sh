#!/bin/bash

# Run this file to produce Fig. 4 a) and b) in the paper.
# Set any of the reruns below to yes in order to rerun and recreate that dataset
# Otherwise, use already computed data stored below in ascii.
# This bash-file generates a python-file with speed-data, runs it and plots 
# fig 4 a) and b). The figures are stored as speedleg2cheba.png and speedleg2chebb.png. 
# If these files exist they will be overwritten.
# As long as the trusted data are stored in this file it is perfectly safe
# to rerun either the fmm or the dct. The new data will only appear in a 
# new generated python file (see filename below).

rerun_fmm="no"
rerun_dct="no"
rerun_direct="no"
rerun_plan="no"

smax=32
threads=1
filename="speed_USEACC_${USE_ACCELERATE}_s${smax}_threads${threads}_rerunfmm_${rerun_fmm}_rerundct_${rerun_dct}_rerundirect_${rerun_direct}_rerunplan_${rerun_plan}_${RANDOM}.py"
echo $filename
echo "import numpy as np" > $filename
echo import matplotlib.pyplot as plt >> $filename
if [ "$rerun_fmm" = "yes" ]; then
for (( m=0; m<2; m++))
do
echo "time_${m} = np.array([" >> $filename
for (( n=8; n<=23; n++ ))
do
    result=$(echo "2^($n)" | bc)
    repeat=$(echo "2^(26-$n)" | bc)
    output=$(l2c -N$result -d0 -r$repeat -s$smax -l$m -v0 -t$threads)
    number=$(echo "$output" | awk '{print $10}')
    echo "[$result, $number]," >> $filename
    sleep 0.1
done
echo "])" >> $filename
done
else 
cat << EOF >> $filename 
time_0 = np.array(
[[256, 2.3330e-06],
[512, 6.2500e-06],
[1024, 1.5333e-05],
[2048, 3.3875e-05],
[4096, 6.8542e-05],
[8192, 1.3613e-04],
[16384, 2.6425e-04],
[32768, 5.1150e-04],
[65536, 1.0126e-03],
[131072, 2.0101e-03],
[262144, 3.9895e-03],
[524288, 8.2407e-03],
[1048576, 1.6530e-02],
[2097152, 3.4615e-02],
[4194304, 6.9468e-02],
[8388608, 1.4112e-01]]
)
time_1 = np.array(
[[256, 2.3750e-06],
[512, 6.5410e-06],
[1024, 1.8042e-05],
[2048, 3.9208e-05],
[4096, 7.8209e-05],
[8192, 1.5300e-04],
[16384, 3.0183e-04],
[32768, 5.8475e-04],
[65536, 1.1701e-03],
[131072, 2.3312e-03],
[262144, 4.6654e-03],
[524288, 9.5404e-03],
[1048576, 1.9628e-02],
[2097152, 4.0945e-02],
[4194304, 8.1088e-02],
[8388608, 1.6353e-01]]
)
EOF
fi

# direct
if [ "$rerun_direct" = "yes" ]; then
echo "direct_osx = np.array([" >> $filename
for (( n=6; n<15; n++ ))
do
    result=$(echo "2^($n)" | bc)
    repeat=$(echo "2^(20-$n)" | bc)
    output=$(l2c -N$result -d6 -r$repeat -v0)
    number=$(echo "$output" | awk '{print $10}')
    echo "[$result, $number]," >> $filename
    sleep 0.1
done
echo "])" >> $filename
else 
cat << EOF >> $filename
direct_osx = np.array([
[64, 2.9100e-07],
[128, 9.1600e-07],
[256, 2.9160e-06],
[512, 9.7500e-06],
[1024, 3.6750e-05],
[2048, 1.4979e-04],
[4096, 5.8996e-04],
[8192, 2.9732e-03],
[16384, 1.4975e-02]
])
EOF
fi

# dct
if [ "$rerun_dct" = "yes" ]; then
echo "dct = np.array([" >> $filename
for (( n=6; n<=23; n++ ))
do
    result=$(echo "2^($n)" | bc)
    repeat=$(echo "2^(26-$n)" | bc)
    output=$(l2c -N$result -d7 -r$repeat)
    number=$(echo "$output" | awk '{print $11}')
    echo "[$result, $number]," >> $filename
    sleep 0.1
done
echo "])" >> $filename
else 
cat << EOF >> $filename
dct = np.array(
[[256, 3.750000e-07],
[512, 7.080000e-07],
[1024, 1.500000e-06],
[2048, 3.208000e-06],
[4096, 7.458000e-06],
[8192, 1.808300e-05],
[16384, 4.175000e-05],
[32768, 1.055420e-04],
[65536, 2.209170e-04],
[131072, 4.607910e-04],
[262144, 1.066375e-03],
[524288, 2.624000e-03],
[1048576, 6.263459e-03],
[2097152, 1.187904e-02],
[4194304, 2.722250e-02],
[8388608, 6.122371e-02]]
)
EOF
fi

if [ "$rerun_init" = "yes" ]; then
for (( m=0; m<2; m++))
do
echo "plan_${m} = np.array([" >> $filename
for (( n=8; n<=23; n++ ))
do
    result=$(echo "2^($n)" | bc)
    repeat=$(echo "2^(26-$n)" | bc)
    output=$(l2c -N$result -d9 -r$repeat -s$smax -l$m -v0)
    number=$(echo "$output" | awk '{print $10}')
    echo "[$result, $number]," >> $filename
    sleep 0.1
done
echo "])" >> $filename
done
else 
cat << EOF >> $filename 
plan_0 = np.array(
[[256, 6.2500e-06],
[512, 1.7083e-05],
[1024, 4.8292e-05],
[2048, 8.8583e-05],
[4096, 1.9017e-04],
[8192, 3.8738e-04],
[16384, 7.8017e-04],
[32768, 1.5895e-03],
[65536, 3.2136e-03],
[131072, 6.5246e-03],
[262144, 1.3081e-02],
[524288, 2.5918e-02],
[1048576, 5.2884e-02],
[2097152, 1.0569e-01],
[4194304, 2.1292e-01],
[8388608, 4.2631e-01]]
)
plan_1 = np.array(
[[256, 4.8750e-06],
[512, 1.4958e-05],
[1024, 3.1541e-05],
[2048, 5.9625e-05],
[4096, 1.2400e-04],
[8192, 2.4288e-04],
[16384, 4.8492e-04],
[32768, 9.9117e-04],
[65536, 2.0725e-03],
[131072, 4.0476e-03],
[262144, 8.1327e-03],
[524288, 1.6275e-02],
[1048576, 3.2728e-02],
[2097152, 6.5751e-02],
[4194304, 1.3110e-01],
[8388608, 2.6341e-01]]
)
EOF
fi 

# Results from other machines
cat << EOF >> $filename
time_saga = np.array([
[256, 6.0000e-06],
[512, 2.0000e-05],
[1024, 4.1000e-05],
[2048, 8.5000e-05],
[4096, 1.6600e-04],
[8192, 3.3300e-04],
[16384, 7.0400e-04],
[32768, 1.4060e-03],
[65536, 2.8360e-03],
[131072, 5.6930e-03],
[262144, 1.1404e-02],
[524288, 2.3209e-02],
[1048576, 4.6943e-02],
[2097152, 9.4297e-02],
[4194304, 1.8823e-01],
[8388608, 3.7965e-01]
])
DC = np.array([
[512, 0.000057],
[2048, 0.00065],
[8192, 0.0057],
[32768, 0.035],
[131072, 0.19],
[524288, 0.89],
[2097152, 4.02],
[8388608, 19.03]
])
DCP = np.array([
[512, 0.00073],
[2048, 0.0068],
[8192, 0.059],
[32768, 0.48],
[131072, 3.04],
[524288, 18.13],
[2097152, 94.52],
[8388608, 507.21]
])

# Make figure 4 in the paper
NN = np.array(2**np.arange(4, 25))
fig = plt.figure(figsize=(6, 4))
ax1 = fig.gca()
ax1.loglog(time_0[:, 0], time_1[:, 1], 'k', linestyle=(0, (1, 10)), label='FMM (n)')
ax1.loglog(time_1[:, 0], time_0[:, 1], 'k', label='FMM (m)')
ax1.loglog(time_0[:, 0], time_saga[:, 1], 'k-.', label='FMM (m) SAGA')
ax1.loglog(DC[:, 0], DC[:, 1], 'ko', label='TDC [13]')
ax1.loglog(dct[:, 0], dct[:, 1], 'gray', label='DCT')
ax1.loglog(direct_osx[:-2, 0], direct_osx[:-2, 1], 'gray', linestyle='dotted', label='Direct')
ax1.loglog(time_0[:, 0], 0.7*time_0[:, 0]*time_0[3, 1]/time_0[3, 0], 'k--',
           DC[:, 0], np.log(DC[:, 0])**2*DC[:, 0]*1e-8, 'k--',
           NN[3:14], 5e-11*NN[3:14]**2, 'k--')

ax1.set_xlabel('N')
ax1.set_ylabel('Time [s]')
ax1.legend()
ax1.text(10**3, 0.000005, '$\\\mathcal{O}(N)$', rotation=25)
ax1.text(10**6, 0.7, '$\\\mathcal{O}(N(\\\log N)^2)$', rotation=25)
ax1.text(1.2*10**5, 1, '$\\\mathcal{O}(N^2)$', rotation=45)
ax1.set_yticks(np.logspace(-6, 6, 7))
axes = ax1.axis()
ax1.set_ylim(5e-7, 1000)
plt.savefig('figure4a.png')

fig2 = plt.figure(figsize=(6, 4))
ax2 = fig2.gca()
ax2.loglog(
           plan_0[:, 0], plan_0[:, 1], 'k:',
           plan_1[:, 0], plan_1[:, 1], 'k-.',
           DCP[:, 0], DCP[:, 1], 'ko',
           time_0[:, 0], 0.7*time_0[:, 0]*time_0[3, 1]/time_0[3, 0], 'k--',
           DCP[:, 0], np.log(DC[:, 0])**2*DC[:, 0]*1e-8, 'k--',
           )
ax2.set_xlabel('N')
ax2.legend(['FMM (m)', 'FMM (n)', 'TDC [13]'], loc='upper left')
ax2.text(10**3, 0.000005, '$\\\mathcal{O}(N)$', rotation=25)
ax2.text(10**6, 0.7, '$\\\mathcal{O}(N(\\\log N)^2)$', rotation=25)
ax2.set_xlim(axes[:2])
ax2.set_yticks(np.logspace(-6, 6, 7))
ax2.set_ylim(5e-7, 1000)
plt.savefig('figure4b.png')
plt.show()
EOF

python $filename