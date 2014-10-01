addpath ./Study1;

%--------Skapar filtrerat brus------------------------------

N = 2^10;
N2 = N+1;
Ts = 1; %length of the measured signal.
fs = N/Ts; %sampling frequency.
T = Ts/N; %sampling length.

x = randn(1,N);


n = linspace(0,N,N); 
wc = 0.2;
[b,a] = butter(10,wc,'low');

y = filter(b,a,x);
t = linspace(-N/2,N/2,N);
% figure(1)
% subplot(211)
% plot(t,x);xlim([-N/2 N/2]);
% subplot(212)
 plot(t,y);xlim([-N/2 N/2]);
%%
%-----------Skapar system---------------------------

theta0 = (wc-1.5); %vad är ett lämpligt värde?
zsq = y.^2;

zhw = zeros(1,N);
zhw(y>0) = y(y>0);

zam = y.*cos(theta0*n);
%plot av sytemsignaler
figure(2)
subplot(221)
plot(t,zsq);xlim([-N/2 N/2]);
subplot(222)
plot(t,zhw);xlim([-N/2 N/2]);
subplot(223)
plot(t,zam);xlim([-N/2 N/2]);
subplot(224)
plot(t,y);xlim([-N/2 N/2]);
%%
%---------Kolla carrier frequency--------------------
ZAM = fft(zam);
figure(3);
subplot(222);
plot(t,ZAM); title('Histogram of sq');

%%
%----------------PSD av utsignaler----------------
%------------Periodogram----------------------
w = linspace(0,1,N);

Rypsq = PeriodFourier(zsq);
Ryphw = PeriodFourier(zhw);
Rypam = PeriodFourier(zam);

figure(4);
subplot(222);
plot(w,Rypsq); title('Raw Periodogram of sq');
xlabel('[\theta]')
subplot(223);
plot(w,Ryphw);title('Raw Periodogram of hw');
xlabel('[\theta]')
subplot(224);
plot(w,Rypam); title('Raw Periodogram of am');
xlabel('[\theta]')
%%
%----------------PerAv-------------------------------
k = linspace(0,1,2^7);

RypAvsq = PerAv(Rypsq,2^7);
RypAvhw = PerAv(Ryphw,2^7);
RypAvam = PerAv(Rypam,2^7);

figure(5);
subplot(222);
plot(k,RypAvsq); title('Averaged Periodogram of sq');
xlabel('[\theta]')
subplot(223);
plot(k,RypAvhw);title('Averaged Periodogram of hw');
xlabel('[\theta]')
subplot(224);
plot(k,RypAvam); title('Averaged Periodogram of am');
xlabel('[\theta]')

%%
%-------------Smoothed periodogram-------------------
w2 = window(@blackmanharris,65);
w3 = window(@triang,65);

rysq = EstimateACF(zsq,'Blett');
ryhw = EstimateACF(zhw,'Blett');
ryam = EstimateACF(zam,'Blett');

RypSsq = windowing2(rysq,65, w2);
RypSsq = windowing2(rysq,65,w3);

RypShw = windowing2(ryhw,65, w2);
RypShw = windowing2(ryhw,65,w3);

RypSam = windowing2(ryam,65, w2);
RypSam = windowing2(ryam,65,w3);

figure(6);
subplot(222);
plot(w,RypSsq); title('Smoothed Periodogram of sq');
xlabel('[\theta]')
subplot(223);
plot(w,RypShw);title('Smoothed Periodogram of hw');
xlabel('[\theta]')
subplot(224);
plot(w,RypSam); title('Smoothed Periodogram of am');
xlabel('[\theta]')
%%
%--------------Histogram, amp. dist.--------------
L = 2^7;
l = linspace(0,1,L);
[fsq,dsq] = hist(zsq,L);
[fhw,dhw] = hist(zhw,L);
[fam,dam] = hist(zam,L);


figure(7);
subplot(222);
plot(l,fsq); title('Histogram of sq');
subplot(223);
plot(l,fhw);title('Histogram of hw');
subplot(224);
plot(l,fam); title('Histogram of am');
