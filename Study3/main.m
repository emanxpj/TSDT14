addpath ../Study1;

N = 2^12; 

x = randn(1,N);


w = linspace(0,1,N);
wc = 0.1;
[b a] = butter(10,2*wc,'low');

y = filter(b,a,x);

%systems

v1 = ones(1,N);
v1(1:2:end) = 0;
v2 =ones(1,N);
v2(2:2:end) = -1;

yalt = y.*v2;
ydec = y.*v1;


%% 
Ry = PeriodFourier(y);
Ralt = PeriodFourier(yalt);
Rdec  = PeriodFourier(ydec);
w2 = window(@blackmanharris,65);

ry = EstimateACF(y,'Blett');
ralt = EstimateACF(yalt,'Blett');
rdec = EstimateACF(ydec,'Blett');

Ry = windowing2(ry,65, w2);

Ralt = windowing2(ralt,65, w2);

Rdec = windowing2(rdec,65, w2);


figure;
subplot(222);
plot(w,Ry); title('Raw Periodogram of y');
xlabel('[\theta]')
subplot(223);
plot(w,Ralt);title('Raw Periodogram of Alt');
xlabel('[\theta]')
subplot(224);
plot(w,Rdec); title('Raw Periodogram of Dec');
xlabel('[\theta]')

%se sida 176 f�r teori
