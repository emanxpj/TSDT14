addpath ../Study1;
clc
clear
close all;

format long;
N = 2^16; 

x = randn(1,N);


w = linspace(0,1,N);
wc = 0.1;
[b, a] = butter(10,2*wc,'low');

y = filter(b,a,x);
t = linspace(0,N-1,N);

%systems

v1 = cos(pi*t);
v2 = cos(pi*t/2).*cos(pi*t/2);
%v2 = zeros(1,N);
%v2(2:2:end) = 1;

yalt = y.*v1;
ydec = y.*v2;
%% 
% 
% ryMy = EstimateACF(y,'Blett');
% ryMy2 = EstimateACF(yalt,'Blett');
% 
% figure;
% subplot(311);
% plot(t,ryMy); title('Raw Periodogram of y');
% xlabel('[\theta]')
% subplot(312);
% plot(t,ryMy2);title('Raw Periodogram of Alt');
% 
% 

%% 
%Theoretical 

RyT = zeros(1,N);
RyT(abs(w) <wc ) = 1;
RyT(abs(w) >1-wc ) = 1;

RaltT = rectpuls((w-1/2)/(2*wc));

RdecT = 1/4*(rectpuls((w-1/2)/(2*wc))) + 1/4*(rectpuls(w/(2*wc))+rectpuls((w-1)/(2*wc)));
 


%% 
Ry = PeriodFourier(y);
Ralt = PeriodFourier(yalt);
Rdec  = PeriodFourier(ydec);

figure;
subplot(311);
plot(w,Ry); title('Raw Periodogram of y');
hold on; plot(w,RyT,'r'); hold off; 
xlabel('[\theta]')
subplot(312);
plot(w,Ralt);title('Raw Periodogram of Alt');
hold on; plot(w,RaltT,'r'); hold off;
xlabel('[\theta]')
subplot(313);
plot(w,Rdec); title('Raw Periodogram of Dec');
hold on; plot(w,RdecT,'r'); hold off;
xlabel('[\theta]')

print -depsc ../Report/PSD3.eps

w2 = window(@blackmanharris,65);

ry = EstimateACF(y,'Blett');
ralt = EstimateACF(yalt,'Blett');
rdec = EstimateACF(ydec,'Blett');

Ry = windowing2(ry,65, w2);

Ralt = windowing2(ralt,65, w2);

Rdec = windowing2(rdec,65, w2);




figure;
subplot(311);
plot(w,Ry); title('Raw Periodogram of y');
hold on; plot(w,RyT,'r'); hold off; ylim([0 1.2]);
xlabel('[\theta]')
subplot(312);
plot(w,Ralt);title('Raw Periodogram of Alt');
hold on; plot(w,RaltT,'r'); hold off;
xlabel('[\theta]')
subplot(313);
plot(w,Rdec); title('Raw Periodogram of Dec');
hold on; plot(w,RdecT,'r'); hold off;
xlabel('[\theta]')

print -depsc ../Report/PSD3smooth.eps
