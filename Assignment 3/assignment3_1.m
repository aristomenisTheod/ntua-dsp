clc; 
clear all;
close all;

f = 2000;
theta_s = pi/2;

theta = linspace(0, pi, 1000);

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
hold off;

figure(2);
for i = 1:length(microphones)
    subplot(2,2,i);
    semilogy(theta, abs(beam_patterns_N(:,i)));
end

%% Functions 

function res = beta_omega (f, theta, N, d, theta_s)
    c = 340;
    omega = 2 * pi * f;
    part1 = sin((N * omega * d * (cos(theta) - cos(theta_s)))/(2 *c));
    part2 = sin((omega * d * (cos(theta) - cos(theta_s)))/(2 *c));
    
    res = part1/(N * part2);
end