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
disp(H2);
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
%b_new
tol = 0.006;
zeros = roots(b_new);
r = abs(z);
z_out = z(r > 1 + tol);
z_in  = z(r < 1 - tol);
z_uc  = z(abs(r - 1) <= tol);

% Step 4: reflect outside zeros inside
z_out_ref = 1 ./ conj(z_out);
disp(abs(zeros));
disp(z_out_ref);
disp(size(z_in));
disp(size(z_uc));
theta = linspace(0, 2*pi, 500);
plot(cos(theta), sin(theta), '--')   % unit circle
hold on
plot(real(z_in), imag(z_in), 'o')
plot(real(z_uc), imag(z_uc), '.')
grid on
axis equal
xlabel('Real')
ylabel('Imag')
title('Zeros of b\_new in Z-Plane')
legend('Unit Circle', 'Zeros')
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