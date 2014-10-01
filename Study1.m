Ts = 1;
N = 2^10;
dt = Ts/N;
fs = 1/dt;
x=randn(N,1).';

t = linspace(0,Ts,N);
n = linspace(0,100,100);
f = linspace(0,Ts,N);
w = linspace(0,1,N);
a = 0.9;

h1 = (1-a).*a.^n;
H1 = (1-a)./(1-a*exp(-1i*2*pi*w));
X = (1/N)*fft(x);
Y = H1.*X;
y = filter(h1,1,x); %ifft(Y,'symmetric');

subplot(231)
plot(t,x); %plot av bruset
subplot(232)
plot(f,abs(X))
subplot(233)
plot(w,abs(H1))
subplot(234)
plot(f,abs(Y))
subplot(235)
plot(t,y)
subplot(236)
plot(n,h1)
%%
Rx = 1;
Ryt1 = Rx*abs((1-a)./(1-a*exp(-1i*2*pi*w))).^2;
ryt1 = zeros(1,N);
ryt1(1) = (1-a)^2/(1-a^2);

figure(2);

subplot(231)
plot(w,Ryt1); title('Theoretical Ry');
subplot(232)
plot(t,ryt1); xlim([-0.1 0.5]); title('Theoretical ry');
subplot(233)
stem(t,ryt1); title('Theoretical ry');
subplot(234);
ryp = EstimateACF(y, t, 'BmanT');
plot(t,ryp); xlim([-0.1 0.5]); title('Estimated ry BmanT');
subplot(235);
plot(t,ryp); xlim([-0.1 0.5]); title('Estimated ry BmanT');

ryp2 = EstimateACF(y, t, 'Blett');
subplot(236);
plot(t,ryp2); xlim([-0.1 0.5]); title('Estimated ry Blett');

Ryp = Periodogram(y);
%%
figure(3);
subplot(221)
plot(w,Ryt1); title('Theoretical Ry');
subplot(222);
ryp2 = periodogram(y); %matlabs own periodogram.
plot(linspace(0,length(ryp2),length(ryp2)),ryp2); title('Matlab Rypred.')
subplot(223);
plot(w,Ryp);title('Predicted Ry with Periodogram');
hold on;
plot(w,Ryt1,'red');
subplot(224)
ryp3 = PeriodFourier(x);
plot(w,ryp3);

%RypAv = PerAv(y,2^10);
%plot(w,RypAv)

