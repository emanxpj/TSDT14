N = 2^10;
Ts = 1; %length of the measured signal.
fs = N/Ts; %sampling frequency.
T = Ts/N; %sampling length.

x = randn(1,N);


n = linspace(0,N,N); 

a = 0.9;
h1 = (1-a).*a.^n;

y = filter(h1, 1,x);

ryMy = EstimateACF(y,'BmanT');
ryMy2 = EstimateACF(y,'Blett');
ryt1 = zeros(1,N);
ryt1(1) = (1-a)^2/(1-a^2);

figure(1); 
subplot(321);
plot(linspace(-length(ryt1)/2,length(ryt1)/2,length(ryt1)),ryt1); 
subplot(322);
%stem(linspace(-19,20,40),ryt1(0:40));
subplot(323);
plot(linspace(-length(ryMy)/2+1,length(ryMy)/2,length(ryMy)),ryMy);
subplot(324);
ryMy(N+1)
stem(linspace(-19,20,40),ryMy(N-19:N+20));
subplot(325);
plot(n,ryMy2);
subplot(326);
stem(linspace(0,20,20),ryMy2(1:20));
%% 
%PDF estimation.

w = linspace(0,1,N);
Rx = 1;
RyMy1 = Periodogram(y);
Ryt1 = Rx*abs((1-a)./(1-a*exp(-1i*2*pi*w))).^2;

figure(2);
subplot(221);
plot(w,Ryt1);
subplot(222);
plot(w,RyMy1);
hold on;
plot(w,Ryt1,'red');




