clc; clear; close all;

%% Multiband specs
% Band 1: 0 <= w/pi <= 0.4,   0.3 <= A <= 0.4
% Band 2: 0.5 <= w/pi <= 0.7, 0.95 <= A <= 1
% Band 3: 0.8 <= w/pi <= 1,   0.05 <= A <= 0.45

%% Midband desired amplitudes and ripples
D1 = (0.3 + 0.4)/2;      % 0.35
D2 = (0.95 + 1)/2;       % 0.975
D3 = (0.05 + 0.45)/2;    % 0.25

d1 = (0.4 - 0.3)/2;      % 0.05
d2 = (1.0 - 0.95)/2;     % 0.025
d3 = (0.45 - 0.05)/2;    % 0.20

%% Parks-McClellan data
% firpm frequency vector is normalized so 1 <-> pi
F = [0 0.4 0.5 0.7 0.8 1.0];
A = [D1 D1 D2 D2 D3 D3];

% Weights inversely proportional to ripple
W = [1/d1 1/d2 1/d3];
W = W / min(W);          % normalize -> [4 8 1]

%% Search for minimum order
N = 3;                   % start from small order
found = false;

while ~found
    % Force even order for Type-I linear phase
    if mod(N,2) ~= 0
        N = N + 1;
    end

    h = firpm(N, F, A, W);

    [H,w] = freqz(h,1,4096);
    mag = abs(H);
    x = w/pi;

    % Evaluate specs in the constrained bands
    idx1 = (x >= 0   & x <= 0.4);
    idx2 = (x >= 0.5 & x <= 0.7);
    idx3 = (x >= 0.8 & x <= 1.0);

    ok1 = all(mag(idx1) >= 0.3)  && all(mag(idx1) <= 0.4);
    ok2 = all(mag(idx2) >= 0.95) && all(mag(idx2) <= 1.0);
    ok3 = all(mag(idx3) >= 0.05) && all(mag(idx3) <= 0.45);

    if ok1 && ok2 && ok3
        found = true;
    else
        N = N + 2;
    end
end

fprintf('Minimum order N = %d\n', N);
fprintf('Filter length   = %d\n', N+1);

%% Ideal desired response for approximation-error plot
wplot = linspace(0, pi, 4096);
xplot = wplot/pi;
Hd = zeros(size(wplot));

Hd(xplot >= 0   & xplot <= 0.4) = D1;
Hd(xplot >= 0.5 & xplot <= 0.7) = D2;
Hd(xplot >= 0.8 & xplot <= 1.0) = D3;

%% Actual response
[H,w] = freqz(h,1,4096);
mag = abs(H);
mag_db = 20*log10(max(mag,1e-8));

%% Approximation error
Err = Hd - mag.';

%% Plots
figure('Name','Parks-McClellan Multiband FIR','NumberTitle','off');

subplot(2,2,1);
stem(0:N, h, 'filled');
xlabel('n');
ylabel('h[n]');
title(['Impulse Response, Order N = ' num2str(N)]);
grid on;

subplot(2,2,2);
plot(xplot, Err, 'LineWidth', 1.2);
xlabel('\omega/\pi');
ylabel('Approx. Error');
title('Approximation Error');
grid on;
xlim([0 1]);

subplot(2,2,3);
plot(x, mag_db, 'LineWidth', 1.2);
xlabel('\omega/\pi');
ylabel('Magnitude (dB)');
title('Magnitude Response (dB)');
grid on;
xlim([0 1]);
ylim([-40 5]);

subplot(2,2,4);
plot(x, mag, 'LineWidth', 1.2);
xlabel('\omega/\pi');
ylabel('|H(e^{j\omega})|');
title('Zoom of Magnitude Response');
grid on;
xlim([0 1]);
ylim([0 1.1]);
hold on;

% Band edges
xline(0.4, '--k');
xline(0.5, '--k');
xline(0.7, '--k');
xline(0.8, '--k');

% Spec bounds
yline(0.3,  '--r');
yline(0.4,  '--r');
yline(0.95, '--b');
yline(1.0,  '--b');
yline(0.05, '--m');
yline(0.45, '--m');

legend('|H(e^{j\omega})|', 'Band edges', '', '', '', ...
       'Band1 bounds', '', 'Band2 bounds', '', 'Band3 bounds', '', ...
       'Location', 'best');