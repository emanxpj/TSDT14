function [Ryp] = Periodogram(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = max(size(x));
Ryp = zeros(1,N);
Xf = zeros(1,N);
for k = 0:N-1
    s = 0;
    for n = 0:N-1
        s = s + x(n+1)*exp(-1i*2*pi*k*n/N);
    end
    Xf(k+1) = s; %The approximate fourier
end
Ryp = (1/N).*abs(Xf).^2;
end

