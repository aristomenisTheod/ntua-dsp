clc;
clear;
close all;

%%%     Part 2       %%%

[music_stereo, Fs] = audioread('music-dsp19.wav');
music = (music_stereo(:,1)+music_stereo(:,2))/2;
music_norma = music/max(music);

win_length = 512;
musicWindowed2 = buffer(music_norma, win_length);
musicWindowedNotNormed = buffer(music, win_length);

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
num_windows = size(musicWindowed2, 2);
u_k = zeros(win_length + L - 1, num_windows, M);
for k = 1:M
    for win = 1:num_windows
        u_k(:, win, k) = conv(filterbank(:,k), musicWindowed2(:,win)); 
    end
end

% Undersampling %
new_win_size = ceil((win_length+L-1)/M);
undersampled_u_k = zeros(new_win_size, num_windows, M); 
for k = 1:M
    for win = 1:num_windows
        undersampled_u_k(:, win, k) = u_k((1:M:win_length+L-1), win, k);
    end
end
        
%% Section 2.2 %%
% Quantization %
run('assignment2_1.m');

% Adaptive number of bits for each frame %
adaptive_quantized = zeros(size(undersampled_u_k,1),size(musicWindowedNotNormed,2),size(undersampled_u_k,3));
Bk = zeros(size(musicWindowedNotNormed,2),1);
for win = 1 : size(undersampled_u_k,3)
    for k = 1 : size(undersampled_u_k,2)
        music_broken = unique(musicWindowedNotNormed(:,k));
        R = length(music_broken);
        Xmin = min(music_broken);
        Xmax = max(music_broken);
        Bk(k,1) = ceil(log2(R/min(Tg(:,k)))-1);
        D = (Xmax-Xmin)/(2.0^Bk(k,1));
        num_quant = 2^Bk(k,1);
        quantoms = zeros(num_quant,1);
        quantoms(1) = Xmin;
        for l = 2 : 2^Bk(k,1)
           quantoms(l,1) = quantoms(l-1,1) + D;
        end
        temp_quantoms = zeros(size(quantoms,1),1);
        for i = 1 : size(undersampled_u_k,1)
           temp_quantoms(:,1) = abs(quantoms(:,1) - undersampled_u_k(i,k,win));
           [min_val, index] = min(temp_quantoms);
           adaptive_quantized(i,k,win) = quantoms(index,1);
        end
    end
end
% figure();
% stairs(adaptive_quantized(:,100));

% 8 bits %

NonAdaptive_quantized = zeros(size(undersampled_u_k,1),size(musicWindowedNotNormed,2),size(undersampled_u_k,3));
Bk = zeros(size(musicWindowedNotNormed,2),1);
for win = 1 : size(undersampled_u_k,3)
    for k = 1 : size(undersampled_u_k,2)
        Xmin = -1;
        Xmax = 1;
        Bk(k,1) = 8;
        D = (Xmax-Xmin)/(2.0^Bk(k,1));
        num_quant = 2^Bk(k,1);
        quantoms = zeros(num_quant,1);
        quantoms(1) = Xmin;
        for l = 2 : 2^Bk(k,1)
           quantoms(l,1) = quantoms(l-1,1) + D;
        end
        temp_quantoms = zeros(size(quantoms,1),1);
        for i = 1 : size(undersampled_u_k,1)
           temp_quantoms(:,1) = abs(quantoms(:,1) - undersampled_u_k(i,k,win));
           [min_val, index] = min(temp_quantoms);
           NonAdaptive_quantized(i,k,win) = quantoms(index,1);
        end
    end
end
% figure();
% stairs(NonAdaptive_quantized(:,100));

%% Section 2.3 %%
% Composition %

% Oversampling %

oversampled_music = zeros(M*size(adaptive_quantized,1),size(adaptive_quantized,2), size(adaptive_quantized,3));
for win = 1 : size(oversampled_music,3)
    for k = 1 : size(oversampled_music,2)
        for i = 1 : size(adaptive_quantized,1)
            oversampled_music(i*M,k,win) = adaptive_quantized(i,k,win);
        end
    end
end

% Filtering %

filterbank = zeros(L, M);

for k = 1:M
    for n = 1:L
        filterbank(n, k) = g_k(n, k, M);
    end
end

filtered = zeros(size(oversampled_music,1)+size(filterbank,1)-1,size(oversampled_music,2), size(oversampled_music,3));
for win = 1 : size(oversampled_music,2)
    for k = 1 : size(oversampled_music,3)
        filtered(:, win, k) = conv(filterbank(:,k), oversampled_music(:,win,k)); 
    end
end

% Final reconstruction %
music_almost = zeros(size(filtered,1),size(filtered,2));
for win = 1 : size(filtered,2)
    for k = 1 : size(filtered,3)
        music_almost(:,win) = music_almost(:,win) + filtered(:, win, k);
    end
end
close all;
figure();
plot(music_almost(:,1));

new_music = zeros((size(music_almost,1)-M+1)*size(music_almost,2),1);
for win = 1 : size(music_almost,2)
    indexing = (size(music_almost,1)-M+1)*(win);
    new_music(indexing:indexing+(size(music_almost,1)-M+1),1) = music_almost(M-1:size(music_almost,1),win);
end

figure();
plot(new_music(:,1));
sound(new_music(:,1));
%% Functions %%

function result = g_k(n,k,M)
%myFun - Description
%
% Syntax: result = g_k(n)
%
% Long description
    result = h_k(2 * M - 1 - n, k, M);
    
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
