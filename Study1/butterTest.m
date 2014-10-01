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
plot(t,ryt1);xlim([-N/2 N/2]); title('Theoretical autocorrelationfunction');
subplot(322);
stem(stemT,ryt1(N/2 -19:N/2+20)); title('Theoretical autocorrelationfunction')
subplot(323);
plot(t,ryMy);xlim([-N/2 N/2]); title('Blackman-Tukey Estimate of ACF');
subplot(324);
stem(stemT,ryMy(N/2-19:N/2+20)); title('Blackman-Tukey Estimate of ACF');
subplot(325);
plot(t,ryMy2);xlim([-N/2 N/2]); title('Bartlett Estimate of ACF');
subplot(326);
stem(stemT,ryMy2(N/2-19:N/2+20)); title('Bartlett Estimate of ACF');



%PSD estimation.

w = linspace(-1/2,1/2,N);
Rx = 1;
RyMy1 = PeriodFourier(y);

RyMy1 = RyMy1([N/2+1:N 1:N/2]);
Ryt1 = zeros(1,N);
Ryt1(abs(w) <a/2 ) = 1;
w = linspace(0,1,N);
w2 = window(@blackmanharris,65);
w3 = window(@triang,65);
RyMy2 = windowing2(ryMy2,65, w2);
RyMy3 = PerAv(RyMy1,2^7);
RyMyTri = windowing2(ryMy2,65,w3);

Ryt1 = Ryt1([N/2+1:N 1:N/2]);

k = linspace(0,1,2^7);

figure(2);
subplot(221);
plot(w,Ryt1);title('Theoretical PSD');
xlabel('[\theta]'); ylim([0 2])
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
plot(k,RyMy3); title('Averaged Periodogram with red theoretical overlay');
hold on;
plot(w,Ryt1,'red');
hold off;
xlabel('[\theta]')

r=linspace(-65/2, 65/2, 65);
figure(3)
subplot(221)
plot(r,w2);title('Blackman-Harris window')
xlabel('k')
subplot(222)
plot(w,RyMy2);title('Smoothed Periodogram with B-H window');
hold on;
plot(w,Ryt1,'red');
hold off;
legend('Estimated','Theoretical')
xlabel('\theta')
subplot(223)
plot(r,w3);title('Triangle window')
xlabel('k')
subplot(224)
plot(w,RyMyTri);title('Smoothed Periodogram with triangle window');
hold on;
plot(w,Ryt1,'red');
hold off;
xlabel('\theta')
legend('Estimated','Theoretical')
print -depsc BSmoothed.eps

figure(4)
plot(w,Ryt1); title('Theoretical PSD');ylim([0 1.2]);
xlabel('\theta')
print -depsc BTheoPSD.eps

figure(5)
subplot(121);
plot(t,ryt1);xlim([-N/2 N/2]); title('Theoretical autocorrelationfunction');
xlabel('k')
subplot(122);
stem(stemT,ryt1(N/2 -19:N/2+20)); title('Theoretical autocorrelationfunction');
xlabel('k')
print -depsc BTheoACF.eps

figure(6)
subplot(221);
plot(t,ryMy);xlim([-N/2 N/2]); title('Blackman-Tukey Estimate of ACF');
xlabel('k')
subplot(222);
stem(stemT,ryMy(N/2-19:N/2+20)); title('Blackman-Tukey Estimate of ACF');
xlabel('k')
subplot(223);
plot(t,ryMy2);xlim([-N/2 N/2]); title('Bartlett Estimate of ACF');
xlabel('k')
subplot(224);
stem(stemT,ryMy2(N/2-19:N/2+20));title('Bartlett Estimate of ACF');
xlabel('k')
print -depsc BEstACF.eps


figure(7)
subplot(121);
plot(w,RyMy1); title('Raw Periodogram');
hold on;
plot(w,Ryt1,'red');
hold off;
legend('Estimated','Theoretical')
xlabel('\theta')
subplot(122);
plot(k,RyMy3); title('Averaged Periodogram');
xlabel('\theta')
hold on;
plot(w,Ryt1,'red');
hold off;
legend('Estimated','Theoretical')
print -depsc BPeriodogram.eps
