function [ Ryp ] = PeriodFourier(x)
%PERIODFOURIER Summary of this function goes here
%   Detailed explanation goes here

N = max(size(x));
%x = x.*hanning(length(x))';
%x = x.*hanning(length(x))';

L = N/2;
n = linspace(-L,L,2*L +1);
k = zeros(1,N);
k(abs(n) < L) = (1-abs(n(abs(n)<L)/L));

x = k.*x;


Xf = ifftshift(fft(x)); %The approximate fourier
Ryp = (1/N)*abs(Xf).^2;



end

