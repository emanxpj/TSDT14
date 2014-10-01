function [ RyS ] = SmoothMat(x,n)
%SMOOTHMAT Summary of this function goes here
%   Detailed explanation goes here

rx = EstimateACF(x,'Blett');
N = length(x);
t = linspace(-N/2,N/2,N);
w= zeros(1,N);

w(abs(t)<n) = 1;

rw = rx.*w;
RyS = PeriodFourier(rw);


end

