Ts = 1;
N = 5000;
dt = Ts/N;
fs = 1/dt;
x=randn(N,1).';

t = linspace(0,Ts,N);
n = linspace(0,100,100);
f = linspace(0,Ts,N);
w = linspace(0,1,N);
a = 0.9;

h1 = (1-a).*a.^n;
H1 = (1-a)./(1-a*exp(-i*2*pi*w));
X = (1/N)*fft(x);
Y = H1.*X;
y = filter(h1,1,x); %ifft(Y,'symmetric');

subplot(231)
plot(t,x); %plot av bruset
subplot(232)
plot(f,abs(X))
subplot(233)
plot(w,H1)
subplot(234)
plot(f,abs(Y))
subplot(235)
plot(t,y)
subplot(236)
plot(n,h1)
%%
Rx = 1;
Ryt1 = Rx*abs((1-a)./(1-a*exp(-i*2*pi*w))).^2;
ryt1 = zeros(1,N);
ryt1(1) = (1-a)^2/(1-a^2);

subplot(231)
plot(w,Ryt1)
subplot(232)
plot(t,ryt1); xlim([-0.1 0.5])
subplot(233)
stem(t,ryt1)

ryp = EstimateACF(y, t, 'BmanT');
subplot(235)
plot(t,ryp); xlim([-0.1 0.5])

ryp2 = EstimateACF(y, t, 'Blett');
subplot(236)
plot(t,ryp2); xlim([-0.1 0.5])

[Ryp RypAv] = Periodogram(y);
RypAv = PerAv(y,1000);
figure(2)
subplot(221)
plot(w,Ryt1)
subplot(223)
plot(w,Ryp)
subplot(224)
plot(w,RypAv)
