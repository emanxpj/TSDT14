addpath ../Study1;

%--------Skapar filtrerat brus------------------------------


% Triangel 

N = 2^11;

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
H(wc>1-wc) =1;
y = ifft(H.*fft(x),'symmetric');
t = linspace(-N/2,N/2,N);
%y = filter(b,a,x);
% figure(1)
% subplot(211)
% plot(t,x);xlim([-N/2 N/2]);
% subplot(212)
 plot(t,y);xlim([-N/2 N/2]);
 
 
 %%
 
 Ry = PeriodFourier(y);
 ry = EstimateACF(y,'Blett');
 figure(8);
 plot(t,ry);
 
 
 
%%
%-----------Skapar system---------------------------

theta0 = (2*pi*wc); %vad är ett lämpligt värde?
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


 Ry = PeriodFourier(y);
 rz = EstimateACF(zsq,'Blett');
 figure(9);
 z3 = abs(fft(rz));
 z3(1:5) = 0;
 plot(w,z3);
 


%%
%---------Kolla carrier frequency--------------------

Y = abs(fft(y));
ZSQ = abs(fft(zsq));
ZHW = abs(fft(zhw));
ZAM = abs(fft(zam));
figure(3);
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

RzPsq = PeriodFourier(zsq);
RzPhw = PeriodFourier(zhw);
RzPam = PeriodFourier(zam);

figure(4);
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
p = 2^3;


k = linspace(0,1,N/p);

RzAvsq = PerAv(zsq,p);
RzAvhw = PerAv(zhw,p);
RzAvam = PerAv(zam,p);

figure(5);
subplot(222);
plot(k,RzAvsq); title('Averaged Periodogram of sq');
xlabel('[\theta]')
subplot(223);
plot(k,RzAvhw);title('Averaged Periodogram of hw');
xlabel('[\theta]')
subplot(224);
plot(k,RzAvam); title('Averaged Periodogram of am');
xlabel('[\theta]')

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

figure(6);
subplot(222);
plot(w,RzSsq); title('Smoothed Periodogram of sq');
xlabel('[\theta]')
subplot(223);
plot(w,RzShw);title('Smoothed Periodogram of hw');
xlabel('[\theta]')
subplot(224);
plot(w,RzSam); title('Smoothed Periodogram of am');
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
plot(dsq,fsq); title('Histogram of sq');
subplot(223);
plot(dhw,fhw);title('Histogram of hw');
subplot(224);
plot(dam,fam); title('Histogram of am');
