clc;
clear;
close all;

% Part 3
% Section 3.1
% http://cvsp.cs.ntua.gr/~nassos/resources/speech_course_2004/OnlineSpeechDemos/speechDemo_2004_Part1.html#2
[speech, Fs] = audioread('speech_utterance.wav');

t = length(speech)/Fs;
tt = linspace(1, t, length(speech));
win_len = 0.02 * Fs;
% TODO: Make it a function
win_overlap = win_len - 1;
ham_win = hamming(win_len);
sigFramed = buffer(speech, win_len, win_overlap);
diagWin = diag(sparse(ham_win));
sigWindowed = diagWin * sigFramed;
energyST = sum(sigWindowed.^2);

delay = fix((win_len - 1)/2);

figure();
plot(tt, speech, tt(1:end-delay), energyST(delay+1:end));

% ZCR

diff = abs(sign([speech; 0]) - sign([0; speech]));
sigFramed2 = buffer(diff, win_len, win_overlap);
diagWin = diag(sparse(ham_win));
sigWindowed2 = diagWin * sigFramed2;
ZCR = sum(sigWindowed2);
ZCR_norm = ZCR/max(ZCR);
figure();
plot(tt, speech, tt(1:end-delay +1), ZCR_norm(delay+1:end));
