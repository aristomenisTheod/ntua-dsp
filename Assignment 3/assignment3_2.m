clc; 
clear;
close all;

%% Part 2: Beamforming application 

%% 2.1: Beamforming in simulated signals

N = 7;
d = 0.08;
theta_s = pi/4;
theta = 3*pi/4;
c = 340;

[clean_signal, Fs] = ...
    audioread('Material\MicArraySimulatedSignals\source.wav');

signal_length = length(clean_signal);

sensors = zeros(signal_length, N);

for i = 1:N
    [sensors(:, i), Fs] = ...
        audioread('Material\MicArraySimulatedSignals\sensor_' ...
            + string(i-1) + '.wav');
end

%% A) Delay-and-sum beamforming

%%% Question 1

omega = linspace(0, 2*pi, signal_length) * Fs;
d_constant = exp(-1i * (N-1) * omega * d * cos(theta_s)/ (c * 2));

d_k = zeros(N, signal_length);
for i = 1:N
    d_k(i, :) = d_constant .* (exp(1i * (i-1) * omega * d * cos(theta_s) / c));
end

filter_H = (1/N) * d_k';

sensors_f = fft(sensors);

y_f = zeros(signal_length, N);
for i = 1:N
    y_f(:, i) = filter_H(:,i) .* sensors_f(:,i);
end

y = ifft(y_f);

y_total = sum(y, 2);

y_total_real = real(y_total);

%%% Question 2

t = linspace(0, signal_length/Fs, signal_length);
%Plots
figure();
plot(t, clean_signal, "-r");
hold on;
plot(t, sensors(:,3), "-g");
plot(t, y_total_real, "-b");
hold off;
title("Signals");
legend(["Clean signal", "Recieced signal from microphone 3", "Delay-and-sum beamformer output"], 'Location', 'southeast');
alpha(.5);

figure();
subplot(2,2,1);
plot(t, clean_signal, "-r");
title("Clean signal");
subplot(2,2,2);
plot(t, sensors(:,3), "-g");
title("Recieced signal from microphone 3");
subplot(2,2,3);
plot(t, y_total_real, "-b");
title("Delay-and-sum beamformer output");

%Spectrograms
figure();
spectrogram(clean_signal, [],[],[],Fs,'yaxis');
title("Clean signal");
figure();
spectrogram(sensors(:,3), [],[],[],Fs, 'yaxis');
title("Recieced signal from microphone 3");
figure();
spectrogram(y_total_real, [],[],[],Fs, 'yaxis');
title("Delay-and-sum beamformer output");

%SNRs
snr_input = snr(sensors(:,3),sensors(:,3)-clean_signal);
snr_output = snr(y_total_real,y_total_real-clean_signal);

%Output file
audiowrite('sim_ds.wav', y_total_real, Fs);

%% B) Single channel Wiener filtering
t_begin = Fs*0.47;
t_end = Fs*0.5;

u_t = sensors(:,3)-clean_signal;
frame_wiener = sensors(t_begin:t_end,3);
noise = u_t(t_begin:t_end,1);

%%% Question 1
L = Fs*0.03;
window = L/3;
noverlap = L/6;
nfft = L+1;

power_spectrum_u = pwelch(noise,L,noverlap,L,Fs,'onesided');
power_spectrum_x = pwelch(frame_wiener,L,noverlap,L,Fs,'onesided');
frequency_responce = 1 - power_spectrum_u./power_spectrum_x;



% linspace(0, length(frequency_responce));
freq = linspace(0,8,length(frequency_responce));
figure();
semilogy(freq,power_spectrum_x);
figure();
semilogy(freq,power_spectrum_u);
figure();
plot(freq,db(frequency_responce));

%%% Question 2
power_spectrum_s = pwelch(clean_signal);
% E = mean((abs(power_spectrum_s - frequence_responce .* power_spectrum_s)).^2);
% 
% nsd = E./power_spectrum_x
nsd = (abs(1-frequency_responce)).^2;

figure();
semilogy(freq,nsd);

%% Beaforming in real signals
N = 7;
d = 0.04;
theta_s = pi/4;

[clean_signal,Fs] = audioread('Material/MicArrayRealSignals/source.wav');
signal_length = length(clean_signal);

sensors = zeros(signal_length, N);

for i = 1:N
    [sensors(:, i), Fs] = ...
        audioread('Material\MicArrayRealSignals\sensor_' ...
            + string(i-1) + '.wav');
end

%% A)Delay and sum beaforming
% Question 1

omega = linspace(0, 2*pi, signal_length) * Fs;
d_constant = exp(-1i * (N-1) * omega * d * cos(theta_s)/ (c * 2));

d_k = zeros(N, signal_length);
for i = 1:N
    d_k(i, :) = d_constant .* (exp(1i * (i-1) * omega * d * cos(theta_s) / c));
end

filter_H = (1/N) * d_k';

sensors_f = fft(sensors);

y_f = zeros(signal_length, N);
for i = 1:N
    y_f(:, i) = filter_H(:,i) .* sensors_f(:,i);
end

y = ifft(y_f);

y_total = sum(y, 2);

y_total_real = real(y_total);

% Question 2
t = linspace(0, signal_length/Fs, signal_length);
%Plots
figure();
plot(t, clean_signal, "-r");
hold on;
plot(t, sensors(:,3), "-g");
plot(t, y_total_real, "-b");
hold off;
title("Signals");
legend(["Clean signal", "Recieced signal from microphone 3", "Delay-and-sum beamformer output"], 'Location', 'southeast');
alpha(.5);

figure();
subplot(2,2,1);
plot(t, clean_signal, "-r");
title("Clean signal");
subplot(2,2,2);
plot(t, sensors(:,3), "-g");
title("Recieced signal from microphone 3");
subplot(2,2,3);
plot(t, y_total_real, "-b");
title("Delay-and-sum beamformer output");

%Spectrograms
figure();
spectrogram(clean_signal, [],[],[],Fs,'yaxis');
title("Clean signal");
figure();
spectrogram(sensors(:,3), [],[],[],Fs, 'yaxis');
title("Recieced signal from microphone 3");
figure();
spectrogram(y_total_real, [],[],[],Fs, 'yaxis');
title("Delay-and-sum beamformer output");

% Question 3



%% FUNCTIONS
