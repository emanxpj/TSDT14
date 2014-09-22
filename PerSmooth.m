function [y2] = PerSmooth(x,N)
%PERSMOOTH Summary of this function goes here
%   Detailed explanation goes here

p = max(size(x));
n = linspace(-p/2,p/2,p+1);
Xf = zeros(1,p);

w1 = zeros(1,p+1);  %y-axis vector initialized to 0, also 500 points like the x-axis vector
w1(abs(n) < N) = 1- abs(n(abs(n)<N)/N); %the points corresponding to |x|< 1 are set to |x|
figure(4);
plot(n,w1);

for k = -p/2:p/2-1
    s = 0;
    for n = 0:p-1
        s = s + w1(k+p/2+1)*x(n+1)*exp(-1i*2*pi*k*n/p);
    end
    Xf(k+1+p/2) = s; %The approximate fourier
end
y2 = 1/(p)*abs(Xf).^2;


end

