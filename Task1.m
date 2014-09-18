%TSDT14 JENS OCH DAVID task1
%% Jens Path
cd ~/Programmering/TSDT14/
%% 1st order lowpassfilter

x=randn(2^10,1).';  
%gaussian noise, if it has constant power it is white. The ACF rx[t1,t2]
%for any t1 != t2, should be zero. rx(t1,t2)  = E[X(t1)X(t2)], these are
%independent with mean zero. For t1 = t2, we get the variance, i.e. 1.
%the power of the noise is thus 1. 

Rx = 1;


%frekvens- och tidsvektor
f = linspace(0,2^10,2^10);
t = linspace(0,2^10,2^10);
dt = 1;
figure(2);
plot(t,x)
%PSD of white gaussian noise is a constant. We use 1 as an example, since
%it leads to simple scalings with different Rx.

%%

%cutoff frequency fc
fc = 10;

%Basic lowpassfilter (1/(1+j*f/fc))
H1 = 1./(1+1i*f/fc);

%Theoretical functions. 
Ryt1 = Rx./(1+(f/fc).^2);
ryt1 = Rx*pi*fc*exp(-2*pi*fc*abs(t));

%Estimated functions, should not use ifft, but a method instead, these are
%simply placeholders.
Ry1 = abs(H1).^2 * Rx;
%ry1 = ifft(Ry1,'symmetric')/dt; %symmetric, due to real signal?
X =fft(x); %X(f);
Y = X.*H1;
y = ifft(Y,'symmetric');
plot(f,y);
%%
ry1 = BmanT(y,t);

figure(1)
subplot(231);
plot(f,Ryt1); xlim([0,10]);
title('Theoretical PSD');
subplot(232);
stem(t,ryt1); xlim([-0.2,20.2]);
title('Theoretical ACF');
subplot(233);
plot(t,ryt1);
title('Theoretical ACF');

subplot(234);
plot(f,Ry1);
title('Estimated PSD');
subplot(235);
stem(t,ry1); xlim([-0.2,20.2]);
title('Estimated ACF');
subplot(236);
plot(t,ry1);
title('Estimated ACF');

%%
%Tio ordningens Butter H

%frekvens- och tidsvektor
f = 0:0.01:99.99;
t = 0:0.01:99.99;
dt = 0.01;
%PSD av vitt Gaussian brus
Rx = 10;
fc = 10;

Wc = 2*pi*fc;
[b,a] = butter(10,2*Wc,'s');

Ry2 = (polyval(b,f)./polyval(a,f))* Rx;
ry2 = ifft(Ry2,'symmetric')/dt;

%Creating the rect
x = linspace(0,99.99,10000);
rect = zeros(size(x));
rect(abs(x)<fc) = 1;


Ryt2 = Rx.*rect;
ryt2 = 2*fc*sinc(2*fc*t)*Rx;

figure(2)
subplot(221);
plot(f,Ryt2);
title('Theoretical PSD');
subplot(222);
plot(t,ryt2); xlim([0,0.5]);
title('Theoretical ACF');

subplot(223);
plot(f,Ry2);
title('Estimated PSD');
subplot(224);
plot(t,ry2); xlim([0,0.5]);
title('Estimated ACF');