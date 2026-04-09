%Sample at Hz and generate a discrete-time signal with
%N=5000 samples in MATLAB (The frequencies, including the sampling frequency, are chosen to make the
%sampling process as arbitrarily as possible). Then quantize the signal by rounding the samples to D=7 bits.
%The rounding can be done by this command in MATLAB: round(x*2^D)/2^D. Find the error between
%the original sampled signal and its quantized version. This error would be quantization error (noise).

%1. Find the power of the quantization error by using the var function in MATLAB. The var
%function calculates the variance of a set of data (or a signal) which is equivalent to the average
%power of the signal. Compare the calculated power of the error with its theoretical value given by
%the fact that the quantization error can be treated as a random variable with uniform distribution
%function. Explain the slight difference between the theoretical value and the calculated one.

%2. Then find the autocorrelation of the error signal. You can use xcorr command in MATLAB to
%do that. Remember that the result should be divided by N because the autocorrelation of a noise
%signal is actually its expected value. Plot the autocorrelation and show that it is an impulse, which
%proves that the error is actually a white noise. The value of the autocorrelation at is the
%average power of the error signal and it should be equal to the value calculated in the previous
%section.

%3. To show that the probability density function of the sample values of error signal is uniform, we
%can use a histogram. You can use hist function in MATLAB to do this. Explain why the
%distribution of the values is not exactly uniform.

%4. Increase the value of N to 30000 and repeat the above steps. Explain the changes in the values.

%5. Increase the number of bits D to 12. Explain the changes in the values.

fs = 29927;

N = 5000;
D = 12; %bits
n = 0:N-1;
t = n / fs;

%% Step 1
x = 0.25*sin(3221*t) + 0.25*cos(8169*t);
xq = round(x * 2^D) / 2^D;
e = x - xq;

step_1_result = var(e);

Delta = 1 / (2^D);
theoretical_power = (Delta^2) / 12;

signal_power = var(x);
SNR_measured_dB = 10*log10(signal_power / step_1_result);
SNR_theoretical_dB = 10*log10(signal_power / theoretical_power);

fprintf('Calculated error power (var): %.10f\n', step_1_result);
fprintf('Theoretical error power:      %.10f\n', theoretical_power);
fprintf('Measured SNR:                 %.10f dB\n', SNR_measured_dB);
fprintf('Theoretical SNR:              %.10f dB\n', SNR_theoretical_dB);
fprintf('Difference:                  %.10f\n', abs(step_1_result - theoretical_power));

figure;
plot(t, x, 'b'); hold on;
plot(t, xq, 'r--');
plot(t, e, 'k');
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal, Quantized Signal, and Error');
legend('x (original)', 'x_q (quantized)', 'e (error)');
grid on;

%% Step 2
[Re, lags] = xcorr(e, 'biased');  % already divides by N
% OR: [Re, lags] = xcorr(e); Re = Re / N;

figure;
stem(lags, Re, 'filled');  % better for impulse-like signals
xlabel('Lag');
ylabel('Autocorrelation');
title('Autocorrelation of Quantization Error');
grid on;

fprintf('Autocorr at lag 0: %.10f\n', Re(lags == 0));
fprintf('Variance of error: %.10f\n', var(e));

%% Step 3
figure;
histogram(e);
grid on;