addpath ../Study1;

close all;
clear all;
clc;

%--------Skapar filtrerat brus------------------------------


% Triangel 

N = 2^12;

Ts = 1; %length of the measured signal.
fs = N/Ts; %sampling frequency.
T = Ts/N; %sampling length.

x = randn(1,N);
w = linspace(0,1,N);

n = linspace(0,N,N); 
wc = 0.1;
[b,a] = butter(10,2*wc,'low');

H = zeros(1,N);
H(w<wc) = 1;
H(w>1-wc) =1;
%y = ifft(H.*fft(x),'symmetric');
t = linspace(-N/2,N/2,N);
y = filter(b,a,x);
% figure(1)
% subplot(211)
% plot(t,x);xlim([-N/2 N/2]);
% subplot(212)
 plot(t,y);xlim([-N/2 N/2]);
 
 
 %%
%.----------------Teoretiska Periodogram---------------

theta0 = (2*pi*(wc+0.2)); %vad är ett lämpligt värde?

Ryt = zeros(1,N);
Ryt(abs(w) <wc ) = 1;
Ryt(abs(w) >1-wc ) = 1;

Rzsqt = 4*wc*(tripuls(w/(4*wc))+tripuls((w-1)/(4*wc)));
Rzsqt(1) = Rzsqt(1) + 4*wc^2;

Rzhwt = 1/(4*pi)*(tripuls(w/(4*wc))+tripuls((w-1)/(4*wc)))+...
    +1/4*(rectpuls(w/(2*wc))+rectpuls((w-1)/(2*wc)));
Rzhwt(1) = Rzhwt(1)+wc/pi;

theta1=theta0/(2*pi);

Rzamt = 1/4*(rectpuls((w-theta1)/(2*wc))+rectpuls((w-1-theta1)/(2*wc))) +...
    + 1/4*(rectpuls((w+theta1)/(2*wc))+rectpuls((w-1+theta1)/(2*wc)));

 
%%
%-----------Skapar system---------------------------

zsq = y.^2;

zhw = zeros(1,N);
zhw(y>0) = y(y>0);

zam = y.*cos(theta0*n);
%plot av sytemsignaler
figure;
subplot(221)
plot(t,zsq);xlim([-N/2 N/2]);
subplot(222)
plot(t,zhw);xlim([-N/2 N/2]);
subplot(223)
plot(t,zam);xlim([-N/2 N/2]);
subplot(224)
plot(t,y);xlim([-N/2 N/2]);

%%


 Ry = PeriodFourier(y);
 rz = EstimateACF(zsq,'Blett');
 figure;
 z3 = abs(fft(rz));
 z3(1:5) = 0;
 plot(w,z3);
 


%%
%---------Kolla carrier frequency--------------------

Y = abs(fft(y));
ZSQ = abs(fft(zsq));
ZHW = abs(fft(zhw));
ZAM = abs(fft(zam));
figure;
subplot(221);
plot(w,Y); title('fft of y');
subplot(222);
plot(w,ZSQ); title('fft of sq');
subplot(223);
plot(w,ZHW); title('fft of hw');
subplot(224);
plot(w,ZAM); title('fft of am');

%%
%----------------PSD av utsignaler----------------
%------------Periodogram----------------------
w = linspace(0,1,N);
Ry = PeriodFourier(y);
RzPsq = PeriodFourier(zsq);
RzPhw = PeriodFourier(zhw);
RzPam = PeriodFourier(zam);

figure;
subplot(222);
plot(w,RzPsq); title('Raw Periodogram of sq');
xlabel('[\theta]')
subplot(223);
plot(w,RzPhw);title('Raw Periodogram of hw');
xlabel('[\theta]')
subplot(224);
plot(w,RzPam); title('Raw Periodogram of am');
xlabel('[\theta]')
%%
%----------------PerAv-------------------------------
% p = 2^3;
% 
% 
% k = linspace(0,1,N/p);
% 
% RzAvsq = PerAv(zsq,p);
% RzAvhw = PerAv(zhw,p);
% RzAvam = PerAv(zam,p);
% 
% figure;
% subplot(222);
% plot(k,RzAvsq); title('Averaged Periodogram of sq');
% xlabel('[\theta]')
% subplot(223);
% plot(k,RzAvhw);title('Averaged Periodogram of hw');
% xlabel('[\theta]')
% subplot(224);
% plot(k,RzAvam); title('Averaged Periodogram of am');
% xlabel('[\theta]')

%%
%-------------Smoothed periodogram-------------------
w2 = window(@blackmanharris,65);
w3 = window(@triang,65);

rzsq = EstimateACF(zsq,'Blett');
rzhw = EstimateACF(zhw,'Blett');
rzam = EstimateACF(zam,'Blett');

RzSsq = windowing2(rzsq,65, w2);
RzSsq2 = windowing2(rzsq,65,w3);

RzShw = windowing2(rzhw,65, w2);
RzShw2 = windowing2(rzhw,65,w3);

RzSam = windowing2(rzam,65, w2);
RzSam2 = windowing2(rzam,65,w3);

figure;
subplot(222);
plot(w,RzSsq); title('Smoothed Periodogram of sq');
hold on; plot(w,Rzsqt,'r');hold off;
xlabel('[\theta]')
subplot(223);
plot(w,RzShw);title('Smoothed Periodogram of hw');
hold on; plot(w,Rzhwt,'r');hold off;
xlabel('[\theta]')
subplot(224);
plot(w,RzSam); title('Smoothed Periodogram of am');
hold on; plot(w,Rzamt,'r');hold off;
xlabel('[\theta]')
%%
%--------------Histogram, amp. dist.--------------
L = 2^4;
l = linspace(0,1,L);
[fsq,dsq] = hist(zsq,L);
[fhw,dhw] = hist(zhw,L);
[fam,dam] = hist(zam,L);

figure;
subplot(221);
hist(y,L); title('Histogram for filtered signal');

subplot(222);
hist(zsq,L);
title('Histogram of sq');
subplot(223);
hist(zhw,L);title('Histogram of hw');
subplot(224);
hist(zam,L); title('Histogram of am');


%%
%Printing of figures.
figure;
subplot(221);
plot(w,Ryt); title('PSD of filtered signal (input to the systems)'); ylim([0 1.25]);
xlabel('\theta');
subplot(222);
plot(w,Rzsqt); title('PSD of halfwave rectified signal'); ylim([0 0.5]);
xlabel('\theta');
subplot(223);
plot(w,Rzhwt); title('PSD of squared signal');ylim([0 0.5]);
xlabel('\theta');
subplot(224);
plot(w,Rzamt); title('PSD of AM-SC modulated signal');ylim([0 0.5]);
xlabel('\theta');

print -depsc ../Report/TheoPSD2.eps

figure;
subplot(221); 
plot(w,Ry); title('Raw Periodogram of input');
xlabel('[\theta]')
hold on; plot(w,Ryt,'r'); hold off;
subplot(222);
plot(w,RzPsq); title('Raw Periodogram of sq');ylim([0 3]);
xlabel('[\theta]')
hold on; plot(w,Rzsqt,'r'); hold off;
subplot(223);
plot(w,RzPhw);title('Raw Periodogram of hw');ylim([0 3]);
xlabel('[\theta]')
hold on; plot(w,Rzhwt,'r'); hold off;
subplot(224);
plot(w,RzPam); title('Raw Periodogram of am');
xlabel('[\theta]')
hold on; plot(w,Rzamt,'r'); hold off;

print -depsc ../Report/RawPSD2.eps

figure;
subplot(221);
hist(y,L); title('Histogram for filtered signal');
subplot(222);
hist(zsq,L);
title('Histogram of sq');
subplot(223);
hist(zhw,L);title('Histogram of hw');
subplot(224);
hist(zam,L); title('Histogram of am');

print -depsc ../Report/Histogram.eps
