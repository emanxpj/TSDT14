N = 2^10;
Ts = 1; %length of the measured signal.
fs = N/Ts; %sampling frequency.
T = Ts/N; %sampling length.

x = randn(1,N);


n = linspace(0,N,N); 

a = 0.8;
h1 = (1-a).*a.^n;

y = filter(h1, 1,x);

ryt1 = (1-a)/(1+a)*a.^(abs(n-(N-1)/2));
 

ryMy = EstimateACF(y,'BmanT');
ryMy2 = EstimateACF(y,'Blett');

stemT = linspace(-19,20,40);
t = linspace(-N/2,N/2,N);
figure(1); 
subplot(321);
plot(t,ryt1);xlim([-N/2 N/2]); title('Theoretical autocorrelationfunction');

subplot(322);
stem(stemT,ryt1(N/2 -19:N/2+20)); title('Theoretical autocorrelationfunction');

subplot(323);
plot(t,ryMy);xlim([-N/2 N/2]); title('Blackman-Tukey Estimate of ACF');

subplot(324);
stem(stemT,ryMy(N/2-19:N/2+20)); title('Blackman-Tukey Estimate of ACF');

subplot(325);
plot(t,ryMy2);xlim([-N/2 N/2]); title('Bartlett Estimate of ACF');

subplot(326);
stem(stemT,ryMy2(N/2-19:N/2+20));title('Bartlett Estimate of ACF');

%PDF estimation.

w = linspace(-1/2,1/2,N);
Rx = 1;
RyMy1 = Periodogram(y);
RyMy1 = RyMy1([N/2+1:N 1:N/2]);
Ryt1 = Rx*abs((1-a)./(1-a*exp(-1i*2*pi*w))).^2;
Ryt1 = Ryt1([N/2+1:N 1:N/2]);
%RyMy2 = PeriodFourier(y);
%RyMy2 = smooth(Periodogram(y),0.1,'loess');
RyMy2 = windowing(y,15,'square');
RyMy2 = RyMy2([N/2+1:N 1:N/2]);

%RyMy2 = PeriodFourier(y);
RyMy3 = PerAv(y,2^7);
RyMy3 = RyMy3([N/2+1:N 1:N/2]);

w = linspace(0,1,N);

figure(2);
subplot(221);
plot(w,Ryt1); title('Theoretical PSD');
xlabel('[\theta]')

subplot(222);
plot(w,RyMy1); title('Raw Periodogram with red theoretical overlay');
hold on;
plot(w,Ryt1,'red');
hold off;
xlabel('[\theta]')

subplot(223);
plot(w,RyMy2);title('Modified (Smoothed) Periodogram with red theoretical overlay');
hold on;
plot(w,Ryt1,'red');
hold off;
xlabel('[\theta]')
subplot(224);
plot(w,RyMy3); title('Averaged Periodogram with red theoretical overlay');
xlabel('[\theta]')
hold on;
plot(w,Ryt1,'red');
hold off;
