addpath ../Study1;

N = 2^12; 

x = randn(1,N);


w = linspace(0,1,N);
wc = 0.1;
[b, a] = butter(10,2*wc,'low');

y = filter(b,a,x);
t = linspace(0,N,N);
%systems

v1 = cos(pi*t);
v2 = cos(pi*t/2).*cos(pi*t/2);

yalt = y.*v1;
ydec = y.*v2;
%% 
%Theoretical 

RyT = zeros(1,N);
RyT(abs(w) <wc ) = 1;
RyT(abs(w) >1-wc ) = 1;

RaltT = 1/2*rectpuls((w-1/2)/(2*wc));

RdecT = 1/8*(rectpuls((w-1/2)/(2*wc))) + 1/4*(rectpuls(w/(2*wc))+rectpuls((w-1)/(2*wc)));
 


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
