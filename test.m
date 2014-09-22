N = 2^10;
Ts = 1; %length of the measured signal.
fs = N/Ts; %sampling frequency.
T = Ts/N; %sampling length.

x = randn(1,N);


n = linspace(0,N,N); 

a = 0.8;
h1 = (1-a).*a.^n;

y = filter(h1, 1,x);

ryMy = EstimateACF(y,'BmanT');
ryMy2 = EstimateACF(y,'Blett');
ryt1 = zeros(1,N);
ryt1(N/2) = (1-a)^2/(1-a^2);

stemT = linspace(-19,20,40);
t = linspace(-N/2,N/2,N);

figure(1); 
subplot(321);
plot(t,ryt1);xlim([-N/2 N/2]); 
subplot(322);
stem(stemT,ryt1(N/2 -19:N/2+20));
subplot(323);
plot(t,ryMy);xlim([-N/2 N/2]);
subplot(324);
stem(stemT,ryMy(N/2-19:N/2+20));
subplot(325);
plot(t,ryMy2);xlim([-N/2 N/2]);
subplot(326);
stem(stemT,ryMy2(N/2-19:N/2+20));

%PDF estimation.

w = linspace(-1/2,1/2,N);
Rx = 1;
RyMy1 = Periodogram(y);
Ryt1 = Rx*abs((1-a)./(1-a*exp(-1i*2*pi*w))).^2;
RyMy2 = PerSmooth2(y,5);
RyMy3 = PerAv(y,2^7);

figure(2);
subplot(221);
plot(w,Ryt1);
subplot(222);
plot(w,RyMy1);
hold on;
plot(w,Ryt1,'red');
hold off;
subplot(223);
plot(linspace(-1/2,1/2,N),RyMy2(1:end));
hold on;
plot(w,Ryt1,'red');
hold off;
subplot(224);
plot(w,RyMy3);
hold on;
plot(w,Ryt1,'red');
hold off;
