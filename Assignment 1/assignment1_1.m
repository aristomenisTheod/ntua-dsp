clear;
clc;
close all;
% Step 1.1

samples = linspace(1,1001,1000);
fs = 8192;

omega_columns = [0.9273, 1.0247, 1.1328];
omega_rows = [0.5346, 0.5906, 0.6535, 0.7217];

tones = cell(10,1);
count = 0;
for row = omega_rows
    for clmn = omega_columns
        
        if (row == 0.7217) && (clmn == 0.9273 || clmn  == 1.1328)
            continue
        else
            count = count + 1;
            tones{count} = createTones(clmn, row, samples);
        end
    end
end

% Step 1.2
D4 = abs(fft(tones{4}));
D6 = abs(fft(tones{6}));

figure();
plot(samples, D4);
title('D4 Plot');

figure();
plot(samples, D6);
title('D6 Plot');

% Step 1.3
% 03113182 + 03115632 = 06228814
pad = zeros(1,100);
AM_sound = [tones{10} pad tones{6} pad tones{2} pad tones{2} pad ...
    tones{8} pad tones{8} pad tones{1} pad tones{4}];
audiowrite('tone_sequence.wav', AM_sound/abs(max(AM_sound)), fs);

figure();
plot(AM_sound);

% Step 1.4
rect_win = rectwin(1000);
hamm_win = hamming(1000);

AM_tones = cell(8, 1);
for i = linspace(1, 8, 8)
    index = 1100*(i-1)+1;
    AM_tones{i} = AM_sound(index: index+999)';
end

F_t_rect = cell(8,1);
F_t_hamm = cell(8,1);

for i = linspace(1, 8, 8)
    F_t_rect{i} = abs(fft(rect_win .* AM_tones{i}));
    F_t_hamm{i} = abs(fft(hamm_win .* AM_tones{i}));
end

for i = linspace(1, 8, 8)
figure();
subplot(1, 2, 1);
plot(F_t_rect{i});
title(sprintf('%s %d', 'Rectangular window digit', i));
subplot(1, 2, 2);
plot(F_t_hamm{i});
title(sprintf('%s %d', 'Hamming window digit', i));
end

% Step 1.5
% https://wiki.analytica.com/index.php?title=FFT
touch_tones = cell(10,2);
touch_tone_freqs = cell(10,2);

for i = linspace(1,10,10)
    tmp = abs(fft(tones{i}));
    tmp = tmp(1:500);
    [~, idxs] = sort(tmp, 'descend');
    peak1 = idxs(1);
    
    for index = idxs
        if abs(index - peak1) < 10
            continue
        else
            peak2 = index;
            break
        end
    end
    touch_tones{i,1} = peak1;
    touch_tones{i,2} = peak2;
    
    touch_tone_freqs{i,1} = peak1*fs/1000;
    touch_tone_freqs{i,2} = peak2*fs/1000;
end
% celldisp(touch_tone_freqs);

% Step 1.6

decoded = ttdecode(AM_sound, touch_tones);

% Step 1.7
load('my_touchtones.mat');
easyDecoded = ttdecode(easySig, touch_tones);
hardDecoded = ttdecode(hardSig, touch_tones); 

function Vector = ttdecode(signIn, touchTones)
    
    [~, c_size] = size(signIn);
    nums = fix(c_size/1000);
    tones = cell(nums,1);
    t_tones = cell(nums, 2);
    hamm_win = hamming(1000)';
    Vector = zeros(1, nums);

    last_index = 0;
    for indx = linspace(1, nums, nums)
        
        index = last_index + 1;
        index_extr = find(signIn(index:end) ~= 0);
        index = index + index_extr(1) - 1;
        last_index = index+999;
        tmp = abs(fft((signIn(index: last_index) .* hamm_win)));
        tones{indx} = tmp(1:500);
        [~, indexes] = sort(tones{indx}, 'descend');
        peak1 = indexes(1);
    
        for inner_index = indexes
            if abs(inner_index - peak1) < 10
                continue
            else
                peak2 = inner_index;
                break
            end
        end
        t_tones{indx, 1} = peak1;
        t_tones{indx, 2} = peak2;
        
        for tt = linspace(1, 10, 10)
            if ((abs(t_tones{indx,1} - touchTones{tt, 1}) < 3) ...
                    && (abs(t_tones{indx,2} - touchTones{tt, 2}) < 3)) 
                if tt == 10
                    Vector(indx) = 0;
                else
                    Vector(indx) = tt;
                end
            end
        end 
    end  
    disp(Vector);
end

function tone = createTones(t1, t2, samples)
    tone = sin(t1*samples) + sin(t2*samples);
end