clear; clc; close all;

% Sample at Hz and generate a discrete-time signal with N=5000 samples
fs = 29927;
N = 5000;
n = 0:N-1;
t = n/fs;

% Input signal
x = 0.1*sin(2*pi*785*t) + 0.1*sin(2*pi*1031*t);

% Distinct cases: [a, D]
cases = [
    0.6   8;
    0.6  12;
    0.95 12
];

for c = 1:size(cases, 1)

    a = cases(c, 1);
    D = cases(c, 2);

    fprintf('\n========================================\n');
    fprintf('Case %d: a = %.2f, D = %d bits\n', c, a, D);
    fprintf('========================================\n');

    % Step 1: Linf scaling
    x_inf = max(abs(x));
    K = (1 - abs(a)) / x_inf;

    fprintf('||x||_inf = %.10f\n', x_inf);
    fprintf('K = %.10f\n', K);

    % Quantizer
    Q = @(v) round(v * 2^D) / 2^D;

    % Quantized input signal
    xq = Q(x);

    % Quantized filter output
    yq = zeros(1, N);

    for i = 2:N
        mult_in = Q(K * xq(i));        % quantized multiplier output for input branch
        mult_fb = Q(a * yq(i-1));      % quantized multiplier output for feedback branch
        yq(i) = Q(mult_in + mult_fb);  % quantized adder output
    end

    figure;
    plot(t, yq, 'b');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Quantized Output of First-Order Digital Filter (a = %.2f, D = %d)', a, D));
    grid on;

    %% Step 4
    y_ideal = zeros(1, N);

    for i = 2:N
        y_ideal(i) = K*x(i) + a*y_ideal(i-1);
    end

    figure;
    plot(t, y_ideal, 'k--', 'LineWidth', 1.2); hold on;
    plot(t, yq, 'b');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Ideal vs Quantized Filter Output (a = %.2f, D = %d)', a, D));
    legend('Ideal output', 'Quantized output');
    grid on;

    e_out = y_ideal - yq;

    figure;
    plot(t, e_out, 'r');
    xlabel('Time (s)');
    ylabel('Error');
    title(sprintf('Quantization Error at Filter Output (a = %.2f, D = %d)', a, D));
    grid on;

    %% Step 5
    % Noise power at output
    noise_power = var(e_out);

    % Ideal output signal power
    signal_power = var(y_ideal);

    % Quantization SNR
    SNR_dB = 10 * log10(signal_power / noise_power);

    fprintf('Ideal output power: %.10f\n', signal_power);
    fprintf('Output noise power: %.10f\n', noise_power);
    fprintf('Quantization SNR:   %.10f dB\n', SNR_dB);

    %% Step 6
    % Quantizer step size
    delta = 2^(-D);
    sigma_q2 = delta^2 / 12;

    % Since this code uses three quantizers inside the loop:
    %   1) Q(K*xq(i))
    %   2) Q(a*yq(i-1))
    %   3) Q(mult_in + mult_fb)
    % and also quantizes the input once:
    %   xq = Q(x)
    %
    % a practical theoretical model is:
    % theoretical_q_power = ((K^2 + 3) * sigma_q2) / (1 - a^2);

    theoretical_q_power = ((K^2 + 3) * sigma_q2) / (1 - a^2);
    theoretical_snr = 10 * log10(signal_power / theoretical_q_power);

    fprintf('Theoretical noise power: %.10f\n', theoretical_q_power);
    fprintf('Theoretical Quantization SNR: %.10f dB\n', theoretical_snr);
    fprintf('Difference (measured - theoretical): %.10f dB\n', SNR_dB - theoretical_snr);

    %% Step 7
    spec_y = fft(y_ideal) / N;   % normalize
    f = (0:N-1) * (fs/N);

    figure;
    plot(f(1:N/2), abs(spec_y(1:N/2)));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(sprintf('Magnitude Spectrum of Output Signal (a = %.2f, D = %d)', a, D));
    grid on;
    xlim([0 fs/2]);

end