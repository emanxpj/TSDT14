function [ PerAv ] = PerAv( x ,k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
N = max(size(x));
l = N/k;
PerAv = zeros(1,max(size(x)));
for p = 1:k
    Xf = zeros(1,l);
    for c = (p-1)*l+1:p*l
        s = 0;
        for n = 1:l
            s =s + x(n)*exp(-1i*2*pi*c*n/N);
        end
        Xf(c) = s;
    end
    PerAv((p-1)*l+1:p*l)=(1/l)*mean(abs(Xf).^2);
end
end

