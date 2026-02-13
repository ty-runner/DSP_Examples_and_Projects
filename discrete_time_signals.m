%VERY USEFUL FUNCTIONS
%1. unit impulse: [x,n] = delta(n1,n0,n2);
%2. unit step: [x,n] = unitstep(n1,n0,n2);
%3. unit pulse: [x,n] = unitpulse(n1,n2,n3,n4);
%4. periodic sequence: x = persegen(xp,Np,Nps)

%HELPERS
%1. Time align, makes sequences same length, [y1,y2,n] =
%timealign(x1,n1,x2,n2);
%2. Folding, reflects a sequence about the origin, [y,ny]=fold(x,nx);
%3. Time shifting, [y,n]=shift(x,n,n0);

n = (-10:10)';

x = 2*cos(2*pi*0.05*n);

disp(x)

figure
stem(n, x, 'fill')
xlabel('n')
ylabel('x[n]')
title('Discrete-Time Signal')
grid on