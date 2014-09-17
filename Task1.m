%TSDT14 JENS OCH DAVID task1

%frekvens- och tidsvektor
f = 0:0.01:99.99;
t = 0:0.01:99.99;
%PSD av vitt Gaussian brus
Rx = 10;
%Simpelt l�gpassfilter H(f) = 1/(1+jf/fc)
%fc �r sk�rfrekvens

H1 = 1./(1+j*f/fc);
Ry1 = abs(H).^2 * Rx;
ry1 = ifft(Ry1);

Ryt1 = Rx./(1+(f/fc).^2);
ryt1 = Rx*2*pi*fc*exp(-2*pi*fc*t);

figure(1)
subplot(221);
title('Estimated PSD');
plot(f,Ryt1);
subplot(222)
title('Estimated ACF');
plot(t,ryt1);

subplot(223);
title('Estimated PSD');
plot(f,Ry1);
subplot(224)
title('Estimated ACF');
plot(t,ry1);

%%
%Tio ordningens Butter H

Wc = 2*pi*fc;

[b,a] = butter(10,Wc,'s');