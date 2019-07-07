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
 plot(t, sensors(:,4), "-g");
 plot(t, y_total_real, "-b");
 hold off;
 title("Signals");
legend(["Clean signal", "Recieced signal from microphone 3",...
    "Delay-and-sum beamformer output"], 'Location', 'southeast');
 alpha(.5);

figure();
subplot(2,2,1);
plot(t, clean_signal, "-r");
title("Clean signal");
subplot(2,2,2);
plot(t, sensors(:,4), "-g");
title("Recieced signal from microphone 3");
subplot(2,2,3);
plot(t, y_total_real, "-b");
title("Delay-and-sum beamformer output");

%Spectrograms
figure();
spectrogram(clean_signal, [],[],[],Fs,'yaxis');
title("Clean signal");
figure();
spectrogram(sensors(:,4), [],[],[],Fs, 'yaxis');
title("Recieced signal from microphone 3");
figure();
spectrogram(y_total_real, [],[],[],Fs, 'yaxis');
title("Delay-and-sum beamformer output");

%SNRs
snr_input = snr(sensors(:,4),sensors(:,4)-clean_signal);
snr_output = snr(y_total_real,y_total_real-clean_signal);

%Output file
audiowrite('sim_ds.wav', y_total_real, Fs);

%% B) Single channel Wiener filtering
t_begin = Fs*0.47;
t_end = Fs*0.5;
 
% Remember: x(t) = u(t) + s(t)
% s -> clean signal
% u -> noise signal

x = sensors(t_begin:t_end,4);
s = clean_signal(t_begin:t_end);
u = x - s;
L = Fs*0.03;
[power_spectrum_u, freq] = pwelch(u, L/3, L/6, L+1, Fs, 'onesided');
[power_spectrum_x, ~] = pwelch(x, L/3, L/6, L+1, Fs, 'onesided');
 
frequency_responce = 1-power_spectrum_u./power_spectrum_x;
figure();
subplot(1,2,1);
hold on;
plot(freq,db(frequency_responce));
title('Wiener Filter Frequency Responce');
xlim([0,8000]);
xlabel('Freq');
ylabel('Amp');
grid on;




%%% Question 2

% E = mean((abs(power_spectrum_s - frequence_responce .* power_spectrum_s)).^2);
% 
% nsd = E./power_spectrum_x
nsd = (abs(1-frequency_responce)).^2;

subplot(1,2,2);
hold on;
plot(freq,db(nsd));
title('Speech Distortion Index');
xlim([0,8000]);
xlabel('Freq'); 
ylabel('Amp');
grid on;
hold off;
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
plot(t, sensors(:,4), "-g");
plot(t, y_total_real, "-b");
hold off;
title("Signals");
legend(["Clean signal", "Recieced signal from microphone 3",...
    "Delay-and-sum beamformer output"], 'Location', 'southeast');
alpha(.5);

figure();
subplot(2,2,1);
plot(t, clean_signal, "-r");
title("Clean signal");
subplot(2,2,2);
plot(t, sensors(:,4), "-g");
title("Recieced signal from microphone 3");
subplot(2,2,3);
plot(t, y_total_real, "-b");
title("Delay-and-sum beamformer output");

%Spectrograms
figure();
spectrogram(clean_signal, [],[],[],Fs,'yaxis');
title("Clean signal");
figure();
spectrogram(sensors(:,4), [],[],[],Fs, 'yaxis');
title("Recieced signal from microphone 3");
figure();
spectrogram(y_total_real, [],[],[],Fs, 'yaxis');
title("Delay-and-sum beamformer output");

% Question 3



%% FUNCTIONS
