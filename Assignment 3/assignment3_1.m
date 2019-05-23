clc; 
clear;
close all;

%% 1.4 Delay-and-sum beam pattern analysis

f = 2000;
theta_s = pi/2;

theta = linspace(0, pi, 1000);

%%% Question 1
microphones = [4, 8, 12, 16];
d = 0.08;

beam_patterns_N = zeros(length(theta), length(microphones));

for mic = 1:length(microphones)
    for t = 1:length(theta)
        beam_patterns_N(t, mic) = ... 
            beta_omega(f, theta(t), microphones(mic), d, theta_s);
    end
end

figure(1);
hold on;
for i = 1:length(microphones)
    semilogy(theta, abs(beam_patterns_N(:,i)));
end
legend(string(microphones));
title('Delay-and-sum beam patterns for different number of microphones');

hold off;

figure(2);
for i = 1:length(microphones)
    subplot(2,2,i);
    semilogy(theta, abs(beam_patterns_N(:,i)));
    title("Delay-and-sum beam patterns for " ...
        + microphones(i) + " microphones");
end

%%% Question 2
N = 8;
d = [0.08, 0.12, 0.16, 0.20];

beam_patterns_d = zeros(length(theta), length(d));

for dist = 1:length(d)
    for t = 1:length(theta)
        beam_patterns_d(t, dist) = ... 
            beta_omega(f, theta(t), N, d(dist), theta_s);
    end
end

figure(3);
hold on;
for i = 1:length(d)
    semilogy(theta, abs(beam_patterns_d(:,i)));
end

legend(string(d));
title('Delay-and-sum beam patterns for different number of distances');

hold off;

figure(4);
for i = 1:length(d)
    subplot(2,2,i);
    semilogy(theta, abs(beam_patterns_d(:,i)));
      title("Delay-and-sum beam patterns for " ...
        + d(i) + " meters distance");
end

%% Functions 

function res = beta_omega (f, theta, N, d, theta_s)
    c = 340;
    omega = 2 * pi * f;
    part1 = sin((N * omega * d * (cos(theta) - cos(theta_s)))/(2 *c));
    part2 = sin((omega * d * (cos(theta) - cos(theta_s)))/(2 *c));
    
    res = part1/(N * part2);
end