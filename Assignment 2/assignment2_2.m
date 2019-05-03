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
        Xmin = min(musicWindowed2(:,k));
        Xmax = max(musicWindowed2(:,k));
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
figure();
stairs(adaptive_quantized(:,100));

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
figure();
stairs(NonAdaptive_quantized(:,100));

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
music_almost(:,1) = music_almost(:,1)/max(music_almost(:,1));
figure();
plot(music_almost(:,1));

new_music = zeros((size(music_almost,1)-2*L+1)*size(music_almost,2),1);
indexing = size(music_almost,1);
new_music(1:size(music_almost,1),1) = ...
        new_music(1:(size(music_almost,1)),1)+music_almost(1:size(music_almost,1),win);
for win = 2 : size(music_almost,2)
    new_music(indexing-2*L:indexing+size(music_almost,1)-2*L-1,1) = ...
        new_music(indexing-2*L:indexing+(size(music_almost,1))-2*L-1,1)+music_almost(1:size(music_almost,1),win);
    indexing = indexing + size(music_almost,1)-L*2;
end
new_music = new_music/max(new_music);
audiowrite('compressed.wav',new_music(1:indexing),44100);
figure();
plot(new_music(1:indexing,1));

% Mean Square Error (MSE) %

% err = immse(music, new_music(1:size(music,1),1));
for i = 1: size(new_music,1)
    if(new_music(i,1) > 0)
        break;
    end
end
err = MSE(music,new_music(i:size(music,1)+i-1,1));
error = sum(err(:,1))/size(err(:,1),1);
figure();
plot(err(:,1));
%% 8bit %%
% Oversampling %

oversampled_music8 = zeros(M*size(NonAdaptive_quantized,1),size(NonAdaptive_quantized,2), size(NonAdaptive_quantized,3));
for win = 1 : size(oversampled_music8,3)
    for k = 1 : size(oversampled_music8,2)
        for i = 1 : size(NonAdaptive_quantized,1)
            oversampled_music8(i*M,k,win) = NonAdaptive_quantized(i,k,win);
        end
    end
end

% Filtering %

filterbank8 = zeros(L, M);

for k = 1:M
    for n = 1:L
        filterbank8(n, k) = g_k(n, k, M);
    end
end

filtered8 = zeros(size(oversampled_music8,1)+size(filterbank8,1)-1,size(oversampled_music8,2), size(oversampled_music8,3));
for win = 1 : size(oversampled_music8,2)
    for k = 1 : size(oversampled_music8,3)
        filtered8(:, win, k) = conv(filterbank8(:,k), oversampled_music8(:,win,k)); 
    end
end

% Final reconstruction %
music_almost8 = zeros(size(filtered8,1),size(filtered8,2));
for win = 1 : size(filtered8,2)
    for k = 1 : size(filtered8,3)
        music_almost8(:,win) = music_almost8(:,win) + filtered8(:, win, k);
    end
end
figure();
plot(music_almost8(:,1));

new_music8 = zeros((size(music_almost8,1)-2*L+1)*size(music_almost8,2),1);
indexing8 = size(music_almost8,1);
new_music8(1:size(music_almost8,1),1) = ...
        new_music8(1:(size(music_almost8,1)),1)+music_almost8(1:size(music_almost8,1),win);
for win = 2 : size(music_almost8,2)
    new_music8(indexing8-2*L:indexing8+size(music_almost8,1)-2*L-1,1) = ...
        new_music8(indexing8-2*L:indexing8+(size(music_almost8,1))-2*L-1,1)+music_almost8(1:size(music_almost8,1),win);
    indexing8 = indexing8 + size(music_almost8,1)-L*2;
end
new_music8 = new_music8/max(new_music8);
audiowrite('compressed8.wav',new_music8(1:indexing8),44100);
figure();
plot(new_music8(1:indexing8,1));

% Mean Square Error (MSE) %

err8 = MSE(music,new_music8(i:size(music,1)+i-1,1));
error8 = sum(err8(:,1))/size(err8(:,1),1);
%% Functions %%
function result = MSE(first, second)
%     result = zeros(size(first,1),1);
    result = (first-second).^2;
end

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
