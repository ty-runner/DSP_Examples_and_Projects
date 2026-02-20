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

[H2,~] = freqz(b,1,1024,fsamp);
Hmag = abs(H2);

% Passband indices
pb = f <= 4000;

% Stopband indices
sb = f >= 4500;

% Actual ripples
delta_p_actual = max(abs(Hmag(pb) - 1));
delta_s_actual = max(Hmag(sb));

fprintf('Actual passband ripple = %f\n', delta_p_actual);
fprintf('Actual stopband ripple = %f\n', delta_s_actual);
%% Plot Both
figure
plot(f,abs(H1),'LineWidth',1.5)
hold on
plot(f,abs(H2),'LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Kaiser','Parks-McClellan')
title('Frequency Response Comparison')