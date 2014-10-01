function [ Ryp ] = PeriodFourier(x)
%PERIODFOURIER Summary of this function goes here
%   Detailed explanation goes here

N = max(size(x));
%x = x.*hanning(length(x))';
%x = x.*hanning(length(x))';
Xf = ifftshift(fft(x)); %The approximate fourier
Ryp = (1/N)*abs(Xf).^2;

Ryp = Ryp([N/2+1:N 1:N/2]);


end

