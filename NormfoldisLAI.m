function [hazn]=NormfoldisLAI(LCT)
delh = 0.0001;
z01  = 0:delh:1;
nz   = length(z01);
%__________________________________________________________________________
% Input Foliage Distribution Parameters: 0 <= zmax <= 1 and 0 < sig < 1
%__________________________________________________________________________
[zmax,sig1,sig0]   =  doub_gaus(LCT);
%__________________________________________________________________________
%
% Compute normalized foliage distribution = haz [h*a(z)]
%__________________________________________________________________________
norm       = zeros(1,nz);
for n      = 1:nz
    if z01(n) <= zmax
        norm(n) = (z01(n)-zmax)/sig0;
    else
        norm(n) = (z01(n)-zmax)/sig1;
    end
end
hazn = exp(-norm.*norm);