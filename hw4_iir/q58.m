wp = 0.45*pi; %0.45 pi radians
Ap = 1; %dB
ws = 0.55*pi; %0.55 pi radians
As = 45; %dB

T = 1;         % sampling period
%Part A: use bilinear transformation and the Chebyshev II approx approach to
%obtain a system function H(z) in the parallel form that satisfies the
%specs
Wp = (2/T)*tan(wp/2);
Ws = (2/T)*tan(ws/2);
[N, Omegac] = cheb2ord(Wp, Ws, Ap, As, 's');
[C,D] = cheby2(N, As, Omegac,'s');

% Bilinear transformation to digital filter
[bz, az] = bilinear(C, D, 1/T);

disp('Digital filter H(z) = B(z)/A(z):');
fprintf('Numerator coefficients bz = \n');
disp(bz);
fprintf('Denominator coefficients az = \n');
disp(az);

%% Parallel form using partial-fraction expansion
[r, p, k] = residuez(bz, az);

disp('Parallel form coefficients from residuez:');
disp('Residues r ='); disp(r);
disp('Poles p ='); disp(p);
disp('Direct term k ='); disp(k);
%Part B: design plots of log magnitude, phase, group delay, and impulse
%responses
nfft = 8192;
[H, w] = freqz(bz, az, nfft);
mag_dB = 20*log10(abs(H) + eps);
disp('mag'); disp(mag_dB);
ph = unwrap(angle(H));
[gd, wg] = grpdelay(bz, az, nfft);

figure;
plot(w/pi, mag_dB, 'LineWidth', 1.5);
grid on;
xlabel('\omega/\pi');
ylabel('Magnitude (dB)');
title('Log Magnitude Response');
hold on;
yline(-Ap, '--r', 'A_p = -1 dB');
yline(-As, '--k', 'A_s = -45 dB');
xline(wp/pi, '--g', '\omega_p');
xline(ws/pi, '--m', '\omega_s');
hold off;

figure;
plot(w/pi, ph, 'LineWidth', 1.5);
grid on;
xlabel('\omega/\pi');
ylabel('Phase (radians)');
title('Phase Response');

figure;
plot(wg/pi, gd, 'LineWidth', 1.5);
grid on;
xlabel('\omega/\pi');
ylabel('Group Delay (samples)');
title('Group Delay');

figure;
impz(bz, az);
grid on;
title('Impulse Response');
%Part C: Determine the exact band edge freqs for the given attentuation
atten = -mag_dB;   % attenuation in dB, positive quantity

% Exact passband edge:
% highest frequency where attenuation is still <= Ap
idx_p = find(atten <= Ap, 1, 'last');

if idx_p < length(w)
    % linear interpolation between idx_p and idx_p+1
    w1 = w(idx_p);     a1 = atten(idx_p);
    w2 = w(idx_p+1);   a2 = atten(idx_p+1);
    wp_exact = w1 + (Ap - a1)*(w2 - w1)/(a2 - a1);
else
    wp_exact = w(idx_p);
end

% Exact stopband edge:
% first frequency where attenuation becomes >= As
idx_s = find(atten >= As, 1, 'first');

if idx_s > 1
    % linear interpolation between idx_s-1 and idx_s
    w1 = w(idx_s-1);   a1 = atten(idx_s-1);
    w2 = w(idx_s);     a2 = atten(idx_s);
    ws_exact = w1 + (As - a1)*(w2 - w1)/(a2 - a1);
else
    ws_exact = w(idx_s);
end

fprintf('Exact passband-edge frequency for %.2f dB attenuation:\n', Ap);
fprintf('wp_exact = %.8f rad/sample = %.8f*pi rad/sample\n', ...
    wp_exact, wp_exact/pi);

fprintf('Exact stopband-edge frequency for %.2f dB attenuation:\n', As);
fprintf('ws_exact = %.8f rad/sample = %.8f*pi rad/sample\n', ...
    ws_exact, ws_exact/pi);
