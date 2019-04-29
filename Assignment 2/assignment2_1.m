clc;
clear;
close all;

% Part 1
% Section 1.0

[music_stereo, Fs] = audioread('music-dsp19.wav');
music = (music_stereo(:,1)+music_stereo(:,2))/2;
music_norm = music/max(music);

t = length(music_norm)/Fs;
tt = linspace(1, t, length(music_norm));
figure();
plot(tt, music_norm);

win_len = 512;
ham_win = hann(win_len);

musicFramed = buffer(music_norm, win_len);
diagWin = diag(sparse(ham_win));
musicWindowed = diagWin * musicFramed;

%tl = linspace(1, 512, 512);
%figure();
%plot(tl, musicWindowed(:,1));

% Section 1.1

musicPK = pk(musicWindowed);
k = linspace(0, 256, 
figure();
plot(musicPK(:,2));

function p = pk(signal)
    sigFft = abs(fft(signal)).^2;
    p = 90.302 + 10*log10(sigFft);
end
function b = bark(f)
    b = 13*atan(0.00076*f)+3.5*atan((f/7500)^2);
end