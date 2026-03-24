clc; clear; close all;

num_lp = 10;
den_lp = [1 1];

wc = 100;
Td = 2;

% Step A
[num_hp, den_hp] = lp2hp(num_lp, den_lp, wc);

z_hp = roots(num_hp);
p_hp = roots(den_hp);

figure;
plot(real(z_hp), imag(z_hp), 'ob', 'LineWidth', 2); hold on;
plot(real(p_hp), imag(p_hp), 'xr', 'LineWidth', 2);
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('Analog Highpass Pole-Zero Plot');
legend('Zeros', 'Poles');

% Step B
[bz, az] = bilinear(num_hp, den_hp, 1/Td);
disp(az)

figure;
zplane(bz, az);
title('Digital Filter Pole-Zero Plot');

% Step C
[H, w] = freqz(bz, az, 2048);
figure;
plot(w, 20*log10(abs(H)+eps), 'LineWidth', 1.5);
grid on;
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Magnitude Response of Digital Highpass Filter');