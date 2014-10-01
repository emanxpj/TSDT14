function [per] = PerAv(y,k)
%x is our filtered noise signal
%Split the function in k segments. 
%calculate the ACF for each segment and the PSD. 
%Average the several PSD's to eachother
%
N = max(size(y));
l = N/k;

for p = 1:k
    PerAv2(1:l,p) = EstimateACF(y((p-1)*l+1:p*l),'Blett');
end

per = abs(fft(mean(PerAv2')));