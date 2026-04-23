clear; clc; close all;

%% Common parameters
fs = 29927;
N  = 5000;
n  = 0:N-1;
t  = n/fs;

x = 0.1*sin(2*pi*785*t) + 0.1*sin(2*pi*1031*t);

% Cases requested in the assignment
% [a, D]
cases = [
    0.6   8
    0.6  12
    0.95 12
];

num_cases = size(cases,1);

% Store results for summary
results = zeros(num_cases, 8);
% columns:
% 1=a, 2=D, 3=SNR_noNS_meas, 4=SNR_NS_meas,
% 5=improvement_meas, 6=SNR_noNS_theory, 7=SNR_NS_theory, 8=improvement_theory

for c = 1:num_cases

    a = cases(c,1);
    D = cases(c,2);

    fprintf('\n====================================================\n');
    fprintf('Case %d: a = %.2f, D = %d bits\n', c, a, D);
    fprintf('====================================================\n');

    %% Step 1: Linf scaling
    x_inf = max(abs(x));
    K = (1 - abs(a)) / x_inf;

    fprintf('||x||_inf = %.10f\n', x_inf);
    fprintf('K         = %.10f\n', K);

    %% Quantizer
    Q = @(v) round(v * 2^D) / 2^D;
    delta = 2^(-D);
    sigma_q2 = delta^2 / 12;

    %% Quantized input
    xq = Q(x);

    %% Ideal output (same reference for both systems)
    y_ideal = zeros(1, N);
    for i = 2:N
        y_ideal(i) = K*x(i) + a*y_ideal(i-1);
    end

    signal_power = var(y_ideal);

    %% =========================================================
    %% A) WITHOUT NOISE SHAPING
    %% =========================================================
    yq_noNS = zeros(1, N);

    % Plain recursive quantized filter:
    % y[n] = Q( Q(K*xq[n]) + Q(a*y[n-1]) )
    for i = 2:N
        mult_in = Q(K * xq(i));
        mult_fb = Q(a * yq_noNS(i-1));
        yq_noNS(i) = Q(mult_in + mult_fb);
    end

    e_out_noNS = y_ideal - yq_noNS;
    noise_power_noNS = var(e_out_noNS);
    SNR_noNS_meas_dB = 10 * log10(signal_power / noise_power_noNS);

    % Theoretical noise power for non-noise-shaped structure
    % input quantization + 3 quantizers inside loop
    theoretical_noise_noNS = ((K^2 + 3) * sigma_q2) / (1 - a^2);
    SNR_noNS_theory_dB = 10 * log10(signal_power / theoretical_noise_noNS);

    %% =========================================================
    %% B) WITH NOISE SHAPING
    %% =========================================================
    yq_NS = zeros(1, N);   % quantized output
    v_NS  = zeros(1, N);   % quantizer input
    e_NS  = zeros(1, N);   % quantizer error

    % Noise-shaping recursion:
    % v[n] = K*xq[n] + a*y[n-1] - e[n-1]
    % y[n] = Q(v[n])
    % e[n] = y[n] - v[n]
    for i = 2:N
        v_NS(i)  = K*xq(i) + a*yq_NS(i-1) - e_NS(i-1);
        yq_NS(i) = Q(v_NS(i));
        e_NS(i)  = yq_NS(i) - v_NS(i);
    end

    e_out_NS = y_ideal - yq_NS;
    noise_power_NS = var(e_out_NS);
    SNR_NS_meas_dB = 10 * log10(signal_power / noise_power_NS);

    % Theoretical noise power for noise-shaped structure
    theoretical_noise_NS = sigma_q2 * ( K^2/(1 - a^2) + 2/(1 + a) );
    SNR_NS_theory_dB = 10 * log10(signal_power / theoretical_noise_NS);

    %% Improvements
    improvement_meas   = SNR_NS_meas_dB   - SNR_noNS_meas_dB;
    improvement_theory = SNR_NS_theory_dB - SNR_noNS_theory_dB;

    %% Print results
    fprintf('\n--- WITHOUT noise shaping ---\n');
    fprintf('Measured noise power     = %.10e\n', noise_power_noNS);
    fprintf('Measured SNR             = %.10f dB\n', SNR_noNS_meas_dB);
    fprintf('Theoretical noise power  = %.10e\n', theoretical_noise_noNS);
    fprintf('Theoretical SNR          = %.10f dB\n', SNR_noNS_theory_dB);

    fprintf('\n--- WITH noise shaping ---\n');
    fprintf('Measured noise power     = %.10e\n', noise_power_NS);
    fprintf('Measured SNR             = %.10f dB\n', SNR_NS_meas_dB);
    fprintf('Theoretical noise power  = %.10e\n', theoretical_noise_NS);
    fprintf('Theoretical SNR          = %.10f dB\n', SNR_NS_theory_dB);

    fprintf('\n--- Improvement due to noise shaping ---\n');
    fprintf('Measured SNR improvement = %.10f dB\n', improvement_meas);
    fprintf('Theory SNR improvement   = %.10f dB\n', improvement_theory);

    %% Store results
    results(c,:) = [a, D, SNR_noNS_meas_dB, SNR_NS_meas_dB, ...
                    improvement_meas, SNR_noNS_theory_dB, ...
                    SNR_NS_theory_dB, improvement_theory];

    %% =========================================================
    %% Plots
    %% =========================================================

    % Quantized outputs comparison
    figure;
    plot(t, y_ideal, 'k--', 'LineWidth', 1.1); hold on;
    plot(t, yq_noNS, 'r');
    plot(t, yq_NS, 'b');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Ideal vs Quantized Outputs (a = %.2f, D = %d)', a, D));
    legend('Ideal', 'Without noise shaping', 'With noise shaping');
    grid on;

    % Error comparison
    figure;
    plot(t, e_out_noNS, 'r'); hold on;
    plot(t, e_out_NS, 'b');
    xlabel('Time (s)');
    ylabel('Error');
    title(sprintf('Output Error Comparison (a = %.2f, D = %d)', a, D));
    legend('Without noise shaping', 'With noise shaping');
    grid on;

    % Spectrum of ideal output
    Y_ideal = fft(y_ideal)/N;
    f = (0:N-1)*(fs/N);

    figure;
    plot(f(1:floor(N/2)), abs(Y_ideal(1:floor(N/2))), 'k');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(sprintf('Magnitude Spectrum of Ideal Output (a = %.2f, D = %d)', a, D));
    grid on;
    xlim([0 fs/2]);

    % Spectrum comparison of quantized outputs
    Y_noNS = fft(yq_noNS)/N;
    Y_NS   = fft(yq_NS)/N;

    figure;
    plot(f(1:floor(N/2)), abs(Y_noNS(1:floor(N/2))), 'r'); hold on;
    plot(f(1:floor(N/2)), abs(Y_NS(1:floor(N/2))), 'b');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(sprintf('Spectrum Comparison (a = %.2f, D = %d)', a, D));
    legend('Without noise shaping', 'With noise shaping');
    grid on;
    xlim([0 fs/2]);

end

%% =========================================================
%% Summary table
%% =========================================================
fprintf('\n\n===================== SUMMARY =====================\n');
fprintf('   a      D    SNR_noNS(meas)   SNR_NS(meas)   Improve(meas)   SNR_noNS(th)   SNR_NS(th)   Improve(th)\n');
fprintf('--------------------------------------------------------------------------------------------------------\n');
for c = 1:num_cases
    fprintf('%5.2f   %2d    %14.4f   %12.4f   %13.4f   %12.4f   %10.4f   %11.4f\n', ...
        results(c,1), results(c,2), results(c,3), results(c,4), ...
        results(c,5), results(c,6), results(c,7), results(c,8));
end

%% Bar plot of measured SNRs
figure;
bar(results(:,3:4));
set(gca, 'XTickLabel', ...
    compose('a=%.2f, D=%d', results(:,1), results(:,2)));
xlabel('Case');
ylabel('SNR (dB)');
title('Measured SNR: Without vs With Noise Shaping');
legend('Without noise shaping', 'With noise shaping');
grid on;

%% Bar plot of measured improvement
figure;
bar(results(:,5));
set(gca, 'XTickLabel', ...
    compose('a=%.2f, D=%d', results(:,1), results(:,2)));
xlabel('Case');
ylabel('SNR Improvement (dB)');
title('Measured SNR Improvement Due to Noise Shaping');
grid on;