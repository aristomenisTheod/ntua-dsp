clc;
clear;
close all;

% Part 1
% Section 1.0

[music, Fs] = audioread('music-dsp19.wav');

music_norm = music/max(music);

t = length(music_norm)/Fs;
tt = linspace(1, t, length(music_norm));
figure();
plot(tt, music_norm);

win_len = 512/Fs;
ham_win = hamming(win_len);

musicFramed = buffer(music_norm, win_len);
diagWin = diag(sparse(ham_win));
musicWindowed = musicFramed * diagWin;

