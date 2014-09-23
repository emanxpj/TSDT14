N = 2^10;
N2 = N+1;
Ts = 1; %length of the measured signal.
fs = N/Ts; %sampling frequency.
T = Ts/N; %sampling length.

x = randn(1,N);


n = linspace(0,N,N); 
wc = 0.4;
[b,a] = butter(10,wc,'low');

y = filter(b,a,x);

a = wc;

ryMy = EstimateACF(y,'BmanT');
ryMy2 = EstimateACF(y,'Blett');
theta0 = a;
ryt1 = theta0*sinc((n-(N-1)/2)*theta0);

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


%%
%PDF estimation.

w = linspace(-1/2,1/2,N);
Rx = 1;
RyMy1 = Periodogram(y);
Ryt1 = zeros(1,N);
Ryt1(abs(w) <a/2 ) = 1;
RyMy2 = PeriodFourier(y);
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
