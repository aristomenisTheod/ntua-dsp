clear;
clc;
close all;

% Step 1.1

fs = 1000;

n = linspace(0, 2, 2*fs);

x = 2*cos(2*pi*70*n) + 3*sin(2*pi*140*n) + 0.15*randn(1,2*fs);

figure();
plot(n,x);

% Step 1.2

[stft, freq, time] = spectrogram(x, 0.04*fs, 0.02*fs, [], fs);
figure();
surf(time, freq, abs(stft));

% Step 1.3

cwts = cwt(x, 'amor');
