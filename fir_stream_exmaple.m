N=20; ni=(0:N-1); x=(3/4).^ni + 0.1*(randn(size(ni)));

M=5; h=ones(1,M)/M;

s = zeros(1,M); %signal memory

for n=1:N
    xin = x(n);
    s(1)=xin;
    yout=h(1)*s(1);
    for m=M:-1:2
        yout=yout+h(m)*s(m);
        s(m)=s(m-1);
    end

    y(n) = yout; % Store the output for the current sample
end

% ---- PLOTTING ----

figure;

subplot(2,1,1);
stem(ni, x, 'filled');
title('Input Signal x[n]');
xlabel('n');
ylabel('x[n]');
grid on;

subplot(2,1,2);
stem(ni, y, 'filled');
title('Filtered Output y[n] (Moving Average)');
xlabel('n');
ylabel('y[n]');
grid on;