function [ Ryp ] = PeriodFourier(x)
%PERIODFOURIER Summary of this function goes here
%   Detailed explanation goes here

N = max(size(x));
Xf = ifftshift(fft(x)); %The approximate fourier
Ryp = (1/N)*abs(Xf).^2;

end

