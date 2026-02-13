% Example vectors
h = [1 2 1];
x = [1 1 1 1];

% Call your convolution function
y = convvec(h, x);

% Time indices
nx = 0:length(x)-1;
nh = 0:length(h)-1;
ny = 0:length(y)-1;

% Plot x
figure;
stem(nx, x, 'filled');
xlabel('n');
ylabel('x[n]');
title('Input Signal x[n]');
grid on;

% Plot h
figure;
stem(nh, h, 'filled');
xlabel('n');
ylabel('h[n]');
title('Impulse Response h[n]');
grid on;

% Plot y
figure;
stem(ny, y, 'filled');
xlabel('n');
ylabel('y[n]');
title('Output y[n] = x[n] * h[n]');
grid on;