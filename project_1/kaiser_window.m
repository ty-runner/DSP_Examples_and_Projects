fsamp = 20000;
fcuts = [4000 4500];
mags = [1 0]; 
devs = [0.1 0.05];

%% Kaiser Design
[n1,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
n1 = n1 + rem(n1,2);
hh = fir1(n1,Wn,ftype,kaiser(n1+1,beta),"noscale");
fprintf('Kaiser Order %d\n', n1);

[H1,f] = freqz(hh,1,1024,fsamp);
Hmag_kaiser = abs(H1);
%% Step 1: Parks-McClellan Design
[n2,fo,ao,w] = firpmord(fcuts,mags,devs,fsamp);
n2 = n2 + rem(n2,2);
b = firpm(n2,fo,ao,w);
fprintf('Parks-McClellan Order %d\n', n2);
fprintf('Parks-McClellan Order %d\n', size(b));
[H2,f] = freqz(b,1,1024,fsamp);
Hmag = abs(H2);

% Passband indices
pb = f <= 4000;

% Stopband indices
sb = f >= 4500;

% Step 2: Actual ripples stop and pass band
delta_p_actual = max(abs(Hmag(pb) - 1));
delta_s_actual = max(Hmag(sb));

% Step 3: Add the stopband delta to h[n] at the center tap (time domain)
M = n2;
center = M/2 + 1;
b_new = b;
b_new(center) = b_new(center) + delta_s_actual;

% Plot original poles and zeros
tol = 1e-8;
z = roots(b_new);
r = abs(z);
figure
zplane(b_new, 1);

% Step 4 and 5: Only keep roots inside and on unit circle, Remove half of the zeros on the unit circle 
on_uc = z(abs(r - 1) <= tol);
inside_uc = z(r < 1 - tol);
n_uc = size(on_uc,1) / 2; %number to remove
z = [on_uc(1:n_uc); inside_uc];

% Step 6: Scale the frequency response with sqrt(a(M/2) / abs(product(z)))
factor = sqrt(b_new(1) / abs(prod(z)));
b_temp = poly(z); % Recreate the coefficients using the poly function
b_temp = b_temp * factor; % Scale as per slides
b_temp = b_temp ./ sqrt(1+delta_s_actual); %final step divide by sqrt(1 + stopband ripple)

% Plot min phase poles and zeros, notice order is now half of original
figure
zplane(b_temp, 1);
[H2_new,f] = freqz(b_temp,1,1024,fsamp);
Hmag_new = abs(H2_new);

% Display change in pass and stop band ripple magnitude
delta_p_kaiser = max(abs(Hmag_kaiser(pb) - 1));
delta_s_kaiser = max(Hmag_kaiser(sb));

fprintf('Actual passband ripple (KAISER)  = %f\n', delta_p_kaiser);
fprintf('Actual stopband ripple (KAISER)  = %f\n', delta_s_kaiser);

fprintf('Actual passband ripple (before) = %f\n', delta_p_actual);
fprintf('Actual stopband ripple (before) = %f\n', delta_s_actual);

delta_p_after = max(abs(Hmag_new(pb) - 1));
delta_s_after = max(Hmag_new(sb));

fprintf('Actual passband ripple (after)  = %f\n', delta_p_after);
fprintf('Actual stopband ripple (after)  = %f\n', delta_s_after);
fprintf('ORDER (after)  = %f\n', size(b_temp));
%% Plot Only Kaiser Windowing
figure
plot(f,abs(H1),'LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Kaiser Windowing')
%% Plot Both
figure
plot(f,abs(H2),'LineWidth',1.5)
hold on
plot(f,abs(H2_new),'LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('H(e^{j\omega T})','H_{min}(e^{j\omega T})')
title('Frequency Response Comparison')

%% ---- Phase + Magnitude Plots for all 3 filters ----
% Recompute (or reuse) responses for consistent plotting
[H1,f]     = freqz(hh,     1, 1024, fsamp);    % Kaiser
[H2,f]     = freqz(b,      1, 1024, fsamp);    % Parks-McClellan (linear phase)
[H2_new,f] = freqz(b_temp, 1, 1024, fsamp);    % Min-phase

% Magnitude in dB (avoid log(0))
eps_db = 1e-12;
H1_mag_db = 20*log10(abs(H1)     + eps_db);
H2_mag_db = 20*log10(abs(H2)     + eps_db);
H3_mag_db = 20*log10(abs(H2_new) + eps_db);

% Unwrapped phase (radians)
H1_phase = unwrap(angle(H1));
H2_phase = unwrap(angle(H2));
H3_phase = unwrap(angle(H2_new));

% ---- Kaiser: Magnitude + Phase ----
figure
subplot(2,1,1)
plot(f, H1_mag_db, 'LineWidth', 1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Kaiser Window: Magnitude Response')

subplot(2,1,2)
plot(f, H1_phase, 'LineWidth', 1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Kaiser Window: Unwrapped Phase')

% ---- Parks-McClellan: Magnitude + Phase ----
figure
subplot(2,1,1)
plot(f, H2_mag_db, 'LineWidth', 1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Parks-McClellan: Magnitude Response')

subplot(2,1,2)
plot(f, H2_phase, 'LineWidth', 1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Parks-McClellan: Unwrapped Phase')

% ---- Min-Phase: Magnitude + Phase ----
figure
subplot(2,1,1)
plot(f, H3_mag_db, 'LineWidth', 1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Minimum Phase: Magnitude Response')

subplot(2,1,2)
plot(f, H3_phase, 'LineWidth', 1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
title('Minimum Phase: Unwrapped Phase')

% ---- Optional: All magnitudes in one plot (dB) ----
figure
plot(f, H1_mag_db, 'LineWidth', 1.5); hold on
plot(f, H2_mag_db, 'LineWidth', 1.5)
plot(f, H3_mag_db, 'LineWidth', 1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Kaiser','Parks-McClellan','Minimum Phase')
title('Magnitude Response Comparison (dB)')

% ---- Optional: All phases in one plot (unwrapped) ----
figure
plot(f, H1_phase, 'LineWidth', 1.5); hold on
plot(f, H2_phase, 'LineWidth', 1.5)
plot(f, H3_phase, 'LineWidth', 1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
legend('Kaiser','Parks-McClellan','Minimum Phase')
title('Phase Response Comparison (Unwrapped)')