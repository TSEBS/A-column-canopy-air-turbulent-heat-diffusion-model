%--------------------------------------------------------------------------
function [ ubc,usb,I1,doh,z0oh,rough,lroug ] = ... 
               calUUSLAII( uz1,uz0,z01,z0h,delh,vkar,usuh,nz,nzet,nexp,cw )
%__________________________________________________________________________
% COMPUTE ---- Composite wind profile [ubc = u(z)/u(h)]
%              Note ubc(z <= z0h) = 0
%__________________________________________________________________________
ubc    = uz1.*uz0;
ubc(1) = 0;
lgl    = z01 >= z0h;
ubc    = lgl.*ubc;
ubx    = max(ubc);
ubc    = ubc/ubx;
%__________________________________________________________________________
% COMPUTE ---- Canopy Reynolds stress profile  
%              Simple empirical function = hyperbolic cosine
%__________________________________________________________________________
B   = 2/(1-exp(-1));
A   = 4.02-B;
adj = A+B*exp(-0.60*nexp);
switch cw
  case 'hypercos'
    usb = cosh(adj*nzet)/cosh(adj*nexp);
  case 'exponent'
    usb = exp(adj*(nzet-nexp));
end
usb = usb/usb(nz);
%__________________________________________________________________________
% COMPUTE ---- d/h and z0/h for whole stand 
%__________________________________________________________________________
% OLD WAY 01/2016 [I1,doh,z0oh,lroug] = caldoh(delh,vkar,usb,usuh);
[I1,doh,z0oh,rough,lroug] = caldohSP(delh,vkar,z01,usb,usuh);
end