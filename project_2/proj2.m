%Plot the frequency response of the designed filter and the above desired function on the same graph.
%You should use “eig” command to find the eigenvector corresponding to the minimum eigenvalue.
%Please note that the outputs of the “eig” command are not necessarily sorted.
% Step 1, eigenfilter from class for type I
% Desired: D(w): 1 from 0 - 4000, nothing from 4000-5000, 0 from 5000-7000,
% nothing from 7000 - 7500, linear increase from 0 to 1 from 7500 - 8500,
% nothing from 8500 - 9000, 0 from 9000 - 12000
clc; clear; close all;

%% Parameters
fs = 24000;
N = 50;
M = N/2;

% Frequency grid (numerical integration)
num_points = 2000;
w = linspace(0, pi, num_points);   % radians

%% Desired response D(w)
f = w * fs / (2*pi);  % convert to Hz

D = zeros(size(w));

for k = 1:length(w)
    if f(k) <= 4000
        D(k) = 1;
    elseif f(k) >= 5000 && f(k) <= 7000
        D(k) = 0;
    elseif f(k) >= 7500 && f(k) <= 8500
        % linear ramp
        D(k) = (f(k) - 7500) / (8500 - 7500);
    elseif f(k) >= 9000
        D(k) = 0;
    else
        % transition regions → don't care much
        D(k) = 0;
    end
end

%% Weighting (IMPORTANT)
W = ones(size(w));

% emphasize important regions
W(f <= 4000) = 10;
W(f >= 5000 & f <= 7000) = 50;
W(f >= 7500 & f <= 8500) = 5;
W(f >= 9000) = 50;

%% Build Q and p
Q = zeros(M+1, M+1);
p = zeros(M+1, 1);

dw = pi / num_points;

for k = 1:length(w)
    % build c(w)
    c = cos((0:M)' * w(k));   % column vector

    Q = Q + W(k) * (c * c') * dw;
    p = p + W(k) * D(k) * c * dw;
end
[V, D] = eig(Q);

% Extract eigenvalues
eigvals = diag(D);
disp(min(eigvals))
% Find index of smallest eigenvalue
[~, idx] = min(eigvals);

% Corresponding eigenvector
a = V(:, idx);
%% Solve for a
a = Q \ p;

%% Recover h(n) from a
h = zeros(N,1);

h(M+1) = a(1);  % center tap

for n = 1:M
    h(M+1+n) = a(n+1)/2;
    h(M+1-n) = a(n+1)/2;
end

%% Frequency response
[H, f_plot] = freqz(h, 1, 2048, fs);

%% Plot
figure;
plot(f_plot, abs(H), 'b', 'LineWidth', 2); hold on;

% Desired overlay
f_des = [0 4000 5000 7000 7500 8500 9000 12000];
d_des = [1 1 0 0 0 1 0 0];
plot(f_des, d_des, 'r--', 'LineWidth', 2);

xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Eigenfilter Design (From Slides)');
legend('Designed Filter', 'Desired');
grid on;
%{ 
Step 2, eigenfilter with a notch added
Design another linear phase type-1 FIR filter with the same frequency response as shown above with
two notches at 4500Hz and 8000Hz. Please note that the first notch is in the first transition band. You
should use the eigenfilter method with linear constraint that was explained in the class. Plot the
magnitude frequency response of the designed filter. To find the null space of a matrix, you can use the
command ”null” in MATLAB.
%}

