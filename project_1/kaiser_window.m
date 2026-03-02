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

%% Parks-McClellan Design
[n2,fo,ao,w] = firpmord(fcuts,mags,devs,fsamp);
n2 = n2 + rem(n2,2);
b = firpm(n2,fo,ao,w);
fprintf('Parks-McClellan Order %d\n', n2);

[H2,f] = freqz(b,1,1024,fsamp);
%disp(H2);
Hmag = abs(H2);

% Passband indices
pb = f <= 4000;

% Stopband indices
sb = f >= 4500;

% Actual ripples
delta_p_actual = max(abs(Hmag(pb) - 1));
delta_s_actual = max(Hmag(sb));

% Add the stopband delta to h[n] at the center tap (time domain)
M = n2;
center = M/2 + 1;
b_new = b;
b_new(center) = b_new(center) + delta_s_actual;

% Plot original zeros vs modified for min phase
tol = 1e-8;
z = roots(b_new);
r = abs(z);
figure;
plot(real(z), imag(z), 'x', 'LineWidth', 2);
hold on;

% Plot unit circle
theta = linspace(0, 2*pi, 500);
plot(cos(theta), sin(theta), 'r--', 'LineWidth', 1.5);

axis equal;
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('Roots before min phase transform');
legend('Roots', 'Unit Circle');
%z(r > 1 + tol) = 1 ./ conj(z(r > 1 + tol));
% Remove half of the zeros on the unit circle, Only keep roots inside and on unit circle

z = z(r - 1 < tol);
disp(z);
z([12 13]) = []; %we can take out either of the pairs on the UC, TODO - cant be manual
r = abs(z);
disp(r);
figure;
plot(real(z), imag(z), 'x', 'LineWidth', 2);
hold on;

% Plot unit circle
theta = linspace(0, 2*pi, 500);
plot(cos(theta), sin(theta), 'r--', 'LineWidth', 1.5);

axis equal;
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('Roots after min phase transform');
legend('Roots', 'Unit Circle');
% Recreate the coefficients using the poly function
b_temp = poly(z);

%gain_correction = sum(b_new) / sum(b_temp);
%disp(gain_correction);
%b_temp = real(gain_correction * b_temp);
disp(sqrt(1 + delta_s_actual));
b_temp = real(b_temp / (sqrt(1 + delta_s_actual)));
[H2_new,f] = freqz(b_temp,1,1024,fsamp);
Hmag_new = abs(H2_new);
fprintf('Actual passband ripple (before) = %f\n', delta_p_actual);
fprintf('Actual stopband ripple (before) = %f\n', delta_s_actual);

delta_p_after = max(abs(Hmag_new(pb) - 1));
delta_s_after = max(Hmag_new(sb));

fprintf('Actual passband ripple (after)  = %f\n', delta_p_after);
fprintf('Actual stopband ripple (after)  = %f\n', delta_s_after);
%% Plot Both
figure
plot(f,abs(H1),'LineWidth',1.5)
hold on
plot(f,abs(H2_new),'LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Kaiser','Parks-McClellan')
title('Frequency Response Comparison')