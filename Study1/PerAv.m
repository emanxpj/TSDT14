function [ PerAv ] = PerAv(x,k)
%x is our filtered noise signal
%Split the function in k segments. 
%calculate the ACF for each segment and the PSD. 
%Average the several PSD's to eachother
%
N = max(size(x));
l = N/k;
PerAv = zeros(l,k);

for p = 1:k
    PerAv(1:l,p) = abs(fft(EstimateACF(x((p-1)*l+1:p*l),'Blett')))';
end

PerAv = mean(PerAv);
