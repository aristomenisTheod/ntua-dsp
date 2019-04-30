clc;
clear;
close all;

%%%     Part 2       %%%

[music_stereo, Fs] = audioread('music-dsp19.wav');
music = (music_stereo(:,1)+music_stereo(:,2))/2;
music_norm = music/max(music);

win_len = 512;
ham_win = hann(win_len);

musicFramed = buffer(music_norm, win_len);
diagWin = diag(sparse(ham_win));
musicWindowed = diagWin * musicFramed;


%%    Section 2.0     %%
% Filterbank %

M = 32;
L = 2 * M;

filterbank = zeros(L, M);

for k = 1:M
    for n = 1:L
        filterbank(n, k) = h_k(n, k, M);
    end
end


%% Section 2.1 %%

% Filtering %
num_windows = size(musicWindowed, 2);
u_k = zeros(win_len + L - 1, num_windows, M);
for k = 1:M
    for win = 1:num_windows
        u_k(:, win, k) = conv(filterbank(:,k), musicWindowed(:,win)); 
    end
end

% Undersampling %
new_win_size = ceil((win_len+L-1)/M);
undersampled_u_k = zeros(new_win_size, num_windows, M); 
for k = 1:M
    for win = 1:num_windows
        undersampled_u_k(:, win, k) = u_k((1:M:win_len+L-1), win, k);
    end
end
        
%% Section 2.2 %%
% Quantization %



%% Functions %%

function result = g_k(n, M)
%myFun - Description
%
% Syntax: result = g_k(n)
%
% Long description
    result = h_k(2 * M - 1 - n);
    
end

function result = h_k(n, k, M)
%h_k - Description
%
% Syntax: result = h_k(n, k)
%s
% Long description
    tmp1 = sin((n + 1/2) * (pi / (2 * M)));
    tmp2 = sqrt(2 / M);
    tmp3 = cos(((2 * n + M + 1) * (2 * k + 1)*pi) / (4 * M));
    result = tmp1 * tmp2 * tmp3;
end
