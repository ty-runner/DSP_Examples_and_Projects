clear; clc; close all;

%% Common parameters
fs = 29927;
N  = 5000;
n  = 0:N-1;
t  = n/fs;

x = 0.1*sin(2*pi*785*t) + 0.1*sin(2*pi*1031*t);

% Cases requested in the assignment
cases = [
    0.6   8
    0.6  12
    0.95 12
];

for c = 1:size(cases,1)

    a = cases(c,1);
    D = cases(c,2);

    fprintf('\n=============================================\n');
    fprintf('Case %d: a = %.2f, D = %d bits\n', c, a, D);
    fprintf('=============================================\n');

    %% Step 1: Linf scaling
    x_inf = max(abs(x));
    K = (1 - abs(a)) / x_inf;

    fprintf('||x||_inf = %.10f\n', x_inf);
    fprintf('K         = %.10f\n', K);

    %% Quantizer
    Q = @(v) round(v * 2^D) / 2^D;
    delta = 2^(-D);
    sigma_q2 = delta^2 / 12;

    %% Step 2 and 3: sample and quantize input, then run noise-shaping system
    xq = Q(x);                % D-bit input signal

    yq = zeros(1, N);         % quantized output
    v  = zeros(1, N);         % quantizer input
    e  = zeros(1, N);         % quantization error of the internal quantizer

    % Noise-shaping recursion:
    % v[n] = K*xq[n] + a*y[n-1] - e[n-1]
    % y[n] = Q(v[n])
    % e[n] = y[n] - v[n]
    for i = 2:N
        v(i)  = K*xq(i) + a*yq(i-1) - e(i-1);
        yq(i) = Q(v(i));
        e(i)  = yq(i) - v(i);
    end

    figure;
    plot(t, yq, 'b');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Quantized Output with Noise Shaping (a = %.2f, D = %d)', a, D));
    grid on;

    %% Step 4: ideal non-quantized system
    y_ideal = zeros(1, N);

    for i = 2:N
        y_ideal(i) = K*x(i) + a*y_ideal(i-1);
    end

    figure;
    plot(t, y_ideal, 'k--', 'LineWidth', 1.1); hold on;
    plot(t, yq, 'b');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Ideal vs Quantized Output (a = %.2f, D = %d)', a, D));
    legend('Ideal output', 'Quantized output');
    grid on;

    %% Step 5: output noise and measured SNR
    e_out = y_ideal - yq;

    figure;
    plot(t, e_out, 'r');
    xlabel('Time (s)');
    ylabel('Error');
    title(sprintf('Output Quantization Error (a = %.2f, D = %d)', a, D));
    grid on;

    noise_power  = var(e_out);
    signal_power = var(y_ideal);

    SNR_meas_dB = 10 * log10(signal_power / noise_power);

    fprintf('Ideal output power      = %.10e\n', signal_power);
    fprintf('Measured noise power    = %.10e\n', noise_power);
    fprintf('Measured output SNR     = %.10f dB\n', SNR_meas_dB);

    %% Step 6: theoretical SNR
    % Total theoretical output noise =
    % input quantization contribution + internal noise-shaped quantizer contribution
    %
    % input quantization contribution:
    % sigma_q2 * K^2 / (1 - a^2)
    %
    % internal shaped quantizer contribution:
    % sigma_q2 * 2 / (1 + a)
    theoretical_noise_power = sigma_q2 * ( K^2/(1 - a^2) + 2/(1 + a) );

    SNR_theory_dB = 10 * log10(signal_power / theoretical_noise_power);

    fprintf('Theoretical noise power = %.10e\n', theoretical_noise_power);
    fprintf('Theoretical output SNR  = %.10f dB\n', SNR_theory_dB);
    fprintf('SNR difference          = %.10f dB\n', abs(SNR_meas_dB - SNR_theory_dB));

    %% Step 7: spectrum of the output signal
    Y = fft(y_ideal)/N;
    f = (0:N-1)*(fs/N);

    figure;
    plot(f(1:floor(N/2)), abs(Y(1:floor(N/2))));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(sprintf('Magnitude Spectrum of Ideal Output (a = %.2f, D = %d)', a, D));
    grid on;
    xlim([0 fs/2]);

    %% Optional: spectrum of quantized output too
    Yq = fft(yq)/N;

    figure;
    plot(f(1:floor(N/2)), abs(Yq(1:floor(N/2))));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(sprintf('Magnitude Spectrum of Quantized Output (a = %.2f, D = %d)', a, D));
    grid on;
    xlim([0 fs/2]);

end