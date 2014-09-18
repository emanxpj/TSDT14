function [ Ryp RypAv ] = Periodogram(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = max(size(x));
s = 0;
Ryp = zeros(1,N);
RypAv = zeros(1,N);
Xf = zeros(1,N);
for k = 1:N
    for n = 1:N
        s =s + x(n)*exp(-i*2*pi*k*n/N);
    end
    Xf(k) = s;
    s = 0;
end
Ryp = (1/N).*abs(Xf).^2;
Xf = fft(x,N);
RypAv = (1/N)*abs(Xf).^2;
end

