clear;
clc;
close all;

% Part 2

% Section 2.1
% Step a

fs = 1000;

n = linspace(0, 2, 2*fs);

x = 2*cos(2*pi*70*n) + 3*sin(2*pi*140*n) + 0.15*randn(1,2*fs);

figure();
plot(n, x);
% Step b

[stft, freq, time] = spectrogram(x, 0.04*fs, 0.02*fs, [], fs);
figure();
surf1 = surf(time, freq, abs(stft));
surf1.EdgeColor = 'none';
colorbar;

% Step c
% https://www.mathworks.com/help/wavelet/ref/cwtft.html

[s, f] = wavescales('morl', fs);
sig = {x, 1/fs};
cwts = cwtft(sig, 'scales', s, 'wavelet', 'morl');
figure();
surf2 = surf(n, f, abs(cwts.cfs));
surf2.EdgeColor = 'none';
colorbar;

% Section 1.2
% Step a

x2 = 1.7*cos(2*pi*90*n) + 0.15*randn(1,2*fs);
x2(0.625 * 1000) = x2(0.625 * 1000) + 3.4;
x2(0.8 * 1000) = x2(0.8* 1000) + 3.4;

figure();
plot(n, x2);

% Step b

[stft, freq, time] = spectrogram(x2, 0.04*fs, 0.02*fs, [], fs);
figure();
contour(time, freq, abs(stft));
colorbar;

% Step c
% https://www.mathworks.com/help/wavelet/ref/cwtft.html

[s, f] = wavescales('morl', fs);
t = linspace(0, 2, 2*fs);
sig = {x2, 1/fs};
cwts = cwtft(sig, 'scales', s, 'wavelet', 'morl');
figure();
contour(t, f, abs(cwts.cfs));
colorbar;