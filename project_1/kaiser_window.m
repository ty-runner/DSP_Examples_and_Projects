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

%% Step 1: Parks-McClellan Design
[n2,fo,ao,w] = firpmord(fcuts,mags,devs,fsamp);
n2 = n2 + rem(n2,2);
b = firpm(n2,fo,ao,w);
fprintf('Parks-McClellan Order %d\n', n2);
fprintf('Parks-McClellan Order %d\n', size(b));
[H2,f] = freqz(b,1,1024,fsamp);
%disp(H2);
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

% Plot original zeros vs modified for min phase
tol = 1e-8;
z = roots(b_new); %plot with x's
r = abs(z);
figure
zplane(b_new, 1);
% Step 4 and 5: Only keep roots inside and on unit circle, Remove half of the zeros on the unit circle 
on_uc = z(abs(r - 1) <= tol);
inside_uc = z(r < 1 - tol);
n_uc = size(on_uc) / 2; %number to remove
z = [on_uc(1:n_uc); inside_uc];
disp(abs(z));

r = abs(z);

% Recreate the coefficients using the poly function
b_temp = poly(z);
figure
zplane(b_temp, 1);
[H2_new,f] = freqz(b_temp,1,1024,fsamp);
Hmag_new = abs(H2_new);
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