fsamp = 20000;
fcuts = [4000 4500];
mags = [1 0]; %passband gain =1, first index, stopband gain =0, second index
devs = [0.1 0.05];

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
n = n + rem(n,2);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),"noscale");
fprintf('Order %d\n', n);
[H,f] = freqz(hh,1,1024,fsamp);
plot(f,abs(H))
grid

