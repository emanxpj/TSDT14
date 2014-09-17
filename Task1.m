%TSDT14 JENS OCH DAVID task1

%frekvens- och tidsvektor
f = 0:0.01:99;
t = 0:0.01:99;
%PSD av vitt Gaussian brus
Rx = 10;
%Simpelt lågpassfilter H(f) = 1/(1+jf/fc)
%fc är skärfrekvens
fc = 10;

H1 = 1./(1+1i*f/fc);
Ry1 = abs(H1).^2 * Rx;
dt = 0.01;
ry1 = ifft(Ry1,'symmetric')/dt; %symmetric, due to real signal?


Ryt1 = Rx./(1+(f/fc).^2);
ryt1 = Rx*pi*fc*exp(-2*pi*fc*abs(t));

figure(1)
subplot(221);
plot(f,Ryt1);
title('Theoreticall PSD');
subplot(222);
plot(t,ryt1); xlim([0,0.1]);
title('Theoreticall ACF');

subplot(223);
plot(f,Ry1);
title('Estimated PSD');
subplot(224);
plot(t,ry1); xlim([0,0.1]);
title('Estimated ACF');

%%
%Tio ordningens Butter H

%frekvens- och tidsvektor
f = 0:0.01:99.99;
t = 0:0.01:99.99;
%PSD av vitt Gaussian brus
Rx = 10;

Wc = 2*pi*fc;
[b,a] = butter(10,Wc,'s');

Ry2 = abs(H2).^2 * Rx;
ry2 = ifftshift(1/(0.01)*ifft(fftshift(Ry2)));

figure(2)
subplot(221);
plot(f,Ryt2);
title('Theoreticall PSD');
subplot(222);
plot(t,ryt2); xlim([0,0.1]);
title('Theoreticall ACF');

subplot(223);
plot(f,Ry2);
title('Estimated PSD');
subplot(224);
plot(t,ry2); xlim([0,0.1]);
title('Estimated ACF');