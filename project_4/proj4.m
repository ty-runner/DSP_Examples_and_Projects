clear; clc; close all;

%% Parameters
fs = 10000;
N = 500;
num_realizations = 50;
Nfft = 1024;

% Filter: H(z) = 1 / (1 - 1.585 z^-1 + 0.96 z^-2)
b = 1;
a = [1 -1.585 0.96];

%% Ideal spectrum
[H, f] = freqz(b, a, Nfft, fs);
ideal_psd = abs(H).^2 / fs;   % ideal PSD in Hz for unit-variance white noise

%% Storage for periodograms
P_all = zeros(Nfft/2+1, num_realizations);

%% Part 2 style: keep one realization for comparison
w1 = randn(1, N);
x1 = filter(b, a, w1);
[P1, f_p] = periodogram(x1, [], Nfft, fs);

%% Repeat for 50 realizations
for k = 1:num_realizations
    w = randn(1, N);          % new white noise realization
    x = filter(b, a, w);      % colored noise output
    [Pxx, f_p] = periodogram(x, [], Nfft, fs);
    P_all(:, k) = Pxx;
end

%% Ensemble average of the periodograms
P_ens_avg = mean(P_all, 2);

%% Plot: one realization vs ensemble average vs ideal
figure;
plot(f_p, 10*log10(P1), 'b'); hold on;
plot(f_p, 10*log10(P_ens_avg), 'm', 'LineWidth', 1.5);
plot(f,   10*log10(ideal_psd), 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('Periodogram: One Realization vs Ensemble Average vs Ideal');
legend('One realization', 'Ensemble average (50)', 'Ideal PSD');
grid on;

%% Optional: linear-scale plot
figure;
plot(f_p, P1, 'b'); hold on;
plot(f_p, P_ens_avg, 'm', 'LineWidth', 1.5);
plot(f, ideal_psd, 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('PSD');
title('Linear Scale: One Realization vs Ensemble Average vs Ideal');
legend('One realization', 'Ensemble average (50)', 'Ideal PSD');
grid on;

%% Optional: estimate variance across realizations at each frequency
var_across_realizations = var(P_all, 0, 2);

figure;
plot(f_p, var_across_realizations, 'k');
xlabel('Frequency (Hz)');
ylabel('Variance');
title('Variance of Periodogram Across 50 Realizations');
grid on;