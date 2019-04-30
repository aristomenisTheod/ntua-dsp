clc;
clear;
close all;

%%%     Part 1       %%%
%%    Section 1.0     %%
% Signal Normalization %

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

tl = linspace(1, 512, 512);
figure();
plot(tl, musicWindowed(:,1));

%%   Section 1.1   %%
% Spectral Analysis %

freq_scale = linspace(20, 20000, 256);
bark_scale = bark(freq_scale);

musicPK = pk(musicWindowed);
k = linspace(1, 256, 256);
figure();
plot(k, musicPK(:,1));

%% Section 1.2 %%
%  Noise & Signal Maskers Identification %

musicSt = zeros(size(musicPK,1), size(musicPK,2));
musicPTM = zeros(size(musicPK,1), size(musicPK,2));
musicPNM = zeros(size(musicPK,1), size(musicPK,2));

for k = 1 : size(musicPK,2)
    musicSt(:,k) = St(musicPK(:,k));
    musicPTM(:,k) = PTM(musicPK(:,k), musicSt(:,k));
    musicPNM(:,k) = findNoiseMaskers(musicPK(:,k), musicPTM(:,k), bark_scale);
end
figure();
plot(musicPNM(:,1));
%% Section 1.3 %%
% Reduction & Reorganization of Maskers %
musicTq = Tq(freq_scale);
figure();
semilogx(freq_scale, musicTq);
grid on;

music_newPTM = zeros(size(musicPK,1), size(musicPK,2));
music_newPNM = zeros(size(musicPK,1), size(musicPK,2));
for k = 1 : size(musicPK,2)
    [music_newPTM(:,k), music_newPNM(:,k)] = checkMaskers(musicPTM(:,k)', musicPNM(:,k)', musicTq, bark_scale); 
end

figure();
plot(music_newPNM(:,1));

%% Section 1.4 %%
% Individual Masking Thresholds %

[music_TTM, music_TNM] = Thresholds(musicPTM, musicPNM, bark_scale);

%% Section 1.5 %%
% Global Masking Threshold %
sumsTtm = sum(10.^0.1*music_TTM);
sumsTnm = sum(10.^0.1*music_TNM);

Tg = zeros(size(musicPTM, 1), size(musicPTM, 2));
for k = 1 : size(Tg, 2)
    for i = 1 : size(Tg, 1)
        Tg(i,k) = 10*log10(10^(0.1*Tq(i))+sumsTtm(1, i, k)+sumsTtm(1, i, k));
    end
end
figure();
plot(Tg(:,1));

%% FUNCTIONS
function [result_Ttm, result_Tnm] = Thresholds(Ptm, Pnm, b)
    result_Ttm = zeros(size(Ptm,1), length(b), size(Ptm, 2));
    result_Tnm = zeros(size(Ptm,1), length(b), size(Ptm, 2));
    SF_Ttm = SF(Ptm, b);
    SF_Tnm = SF(Pnm, b);
    for k = 1 : size(Ptm, 2)
        for j = 1 : size(Ptm,1)
            for i = 1 : size(Ptm,1)
                if( Ptm(j) > 0)
                    result_Ttm(i,j,k) = Ptm(j,k)-0.275*b(j)+SF_Ttm(i,j)-6.025;
                end
                if(Pnm(j) > 0)
                    result_Tnm(i,j,k) = Pnm(j,k)-0.175*b(j)+SF_Tnm(i,j)-2.025;
                end
            end
        end
    end
end
function result = SF(P, b)
    result = zeros(length(b), length(P));
    for j = 1 : length(b)
       if(P(j) > 0.0)
          for i = 1 : length(b)
             if((b(i) >= b(j)-3) && (b(i) <= b(j)+8))
                D = b(i)-b(j);
                if(D < -1)
                    result(i,j) = 17*D-0.4*P(j)+11;
                elseif(D < 0)
                    result(i,j) = (0.4*P(j)+6)*D;
                elseif(D < 1)
                    result(i,j) = -17*D;
                else
                    result(i,j) = (0.15*P(j)-17)*D-0.15*P(j);
                end
             end
          end
       end
    end
end
function result = Tq(f)
    result = 3.64*(f/1000).^(-0.8)-6.5*exp(-0.6*(f/1000-3.3).^2)+10^(-3)*(f/1000).^4;
end
function result = PTM(signal, signalSt)
    result = zeros(length(signal),1);    
    for k = 1 : length(signal)
        if(signalSt(k) == 0)
            result(k) = 0;
        else
            result(k) = 10*log10(10^(0.1*signal(k-1))...
                +10^(0.1*signal(k))+10^(0.1*signal(k+1)));
        end
    end
end
function result = St(signal)
    result = zeros(length(signal),1);
    for k = 1 : length(signal)
        if(k<3 || k>250)
            result(k) = 0; 
        else
            if(signal(k) > signal(k-1) && signal(k) > signal(k+1))
                if(k<63)
                    if((signal(k) > signal(k-2)+7)...
                        && (signal(k) > signal(k+2)+7))
                        result(k) = 1;
                    else
                        result(k) = 0;
                    end
                elseif(k<127)
                    if((signal(k) > signal(k-2)+7)...
                        && (signal(k) > signal(k-3)+7)...
                        && (signal(k) > signal(k+3)+7)...
                        && (signal(k) > signal(k+2)+7))
                        result(k) = 1;
                    else
                        result(k) = 0;
                    end
                else
                    if((signal(k) > signal(k-2)+7)...
                        && (signal(k) > signal(k-3)+7)...
                        && (signal(k) > signal(k-4)+7)...
                        && (signal(k) > signal(k-5)+7)...
                        && (signal(k) > signal(k-6)+7)...
                        && (signal(k) > signal(k+6)+7)...
                        && (signal(k) > signal(k+5)+7)...
                        && (signal(k) > signal(k+4)+7)...
                        && (signal(k) > signal(k+3)+7)...
                        && (signal(k) > signal(k+2)+7))
                        result(k) = 1;
                    else
                        result(k) = 0;
                    end
                end
            else
                result(k) = 0;
            end
        end
    end
end
function p = pk(signal)
    sigFft = abs(fft(signal)).^2;
    temp = 90.302 + 10*log10(sigFft);
    p = temp(1:256,:);
end
function b = bark(f)
    b = 13*atan(0.00076*f)+3.5*atan(f.^2/7500^2);
end