# microphone-Array-with-MATLAB
ADSP lesson project


This issue is considered one of the issues of null steering. We want to place null in the angles 37 and 45, while the desired source enters the array in the direction of the angle 0. This is done with the help of the beamformer h(f). Considering that the number of null angles is smaller than the number of antennas, we can determine the beamformer in such a way that the phase difference caused by the delay of the signal received from different antennas is compensated by the beamformer and at the same time the interference signals from other sources are prevented. We assume that we have N sources, with N < M, that affect the array in directions θ_1≠θ_2≠⋯≠θ_N≠θ_"d" . These sources are considered interference that we want to cancel completely. Mathematically speaking, we can write the constraints equations to eliminate these interferences.
