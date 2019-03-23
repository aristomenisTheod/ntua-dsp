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
energyST = hamm_ste(speech, win_len);

delay = fix((win_len - 1)/2);
figure();
plot(tt, speech, tt(1:end-delay), energyST(delay+1:end));

% ZCR

ZCR_norm = hamm_zcr(speech, win_len);
figure();
plot(tt, speech, tt(1:end-delay +1), ZCR_norm(delay+1:end));

% Section 3.2
[music, Fs] = audioread('music_cut.wav');

t = length(speech)/Fs;
tt = linspace(1, t, length(music));
win_len = 0.02 * Fs;
energyST = hamm_ste(music, win_len);

delay = fix((win_len - 1)/2);
figure();
plot(tt, music, tt(1:end-delay), energyST(delay+1:end));

% ZCR

ZCR_norm = hamm_zcr(music, win_len);
figure();
plot(tt, music, tt(1:end-delay+1), ZCR_norm(delay+1:end));


function energyST_norm = hamm_ste (signal, win_len)
    win_overlap = win_len - 1;
    ham_win = hamming(win_len);
    sigFramed = buffer(signal, win_len, win_overlap);
    diagWin = diag(sparse(ham_win));
    sigWindowed = diagWin * sigFramed;
    energyST = sum(sigWindowed.^2);
    energyST_norm = energyST/max(energyST);
end

function ZCR_norm = hamm_zcr (signal, win_len)
    win_overlap = win_len - 1;
    ham_win = hamming(win_len);
    difference = abs(sign([signal; 0]) - sign([0; signal]));
    sigFramed2 = buffer(difference, win_len, win_overlap);
    diagWin = diag(sparse(ham_win));
    sigWindowed2 = diagWin * sigFramed2;
    ZCR = sum(sigWindowed2);
    ZCR_norm = ZCR/max(ZCR);
end
