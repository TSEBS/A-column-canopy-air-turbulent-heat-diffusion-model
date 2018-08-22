function [c_u,c_pt]=CuCpt(lmo,zm1,zh1,zm2,zh2)
% Note: zm1和zh1可以是一列数，也可以是一个值，
% zm1=1;zh1=1;zm2=zm;zh2=zh;
%  PURPOSE:
%     Compute C_u and C_pt for both unstable and stable surface layer 
%     based on Monin-Obukhov similarity theory. The univeral profile form 
%     was proposed Hogstrom (1996) and the anlytical solution was 
%     developed by Yang (2000)
%     c_u: Frictional velocity /wind speed
%     c_pt: Nondimensional temperature scale
gammam=19.0;
gammah=11.6;
betam=5.3;
betah=8.0;
kv = 0.4;  % Von Karman constant
prantl01 = 1.0;
prantl02 = 0.95; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(lmo);
c_u(1:m,1:n)=NaN;
c_pt(1:m,1:n)=NaN;

if length(zm1)==1   %只是简单的两层高度时的算法
%#######################################################################
%     Unstable case
%#######################################################################
          lmo1=lmo(lmo<0.0);
          xx2   = sqrt(sqrt(1-gammam* zm2./lmo1));   
          xx1   = sqrt(sqrt(1-gammam* zm1./lmo1));
          psim  = 2 * log((1+xx2)./(1+xx1))+ log((1+xx2.*xx2)./(1+xx1.*xx1))- 2*atan(xx2) + 2 * atan(xx1);
          yy2   = sqrt(1-gammah* zh2./lmo1);  
          yy1   = sqrt(1-gammah* zh1./lmo1);
          psih  = 2 * log((1+yy2)./(1+yy1)); 
          uprf  = max(log(zm2/zm1)- psim,0.50*log(zm2/zm1));
          ptprf = max(log(zh2/zh1) - psih,0.33*log(zm2/zm1)); 
          c_u(lmo<0.0)   = kv./ uprf;
          c_pt(lmo<0.0)  = kv./ (ptprf * prantl02 );           
%#######################################################################
%     Stable case:
%#######################################################################
          lmo1=lmo(lmo>=0.0);
          psim  = -betam * (zm2-zm1)./lmo1;
          psih  = -betah * (zh2-zh1)./lmo1;
          psim  = max( -betam, psim);   %实数的对比，fortran用amax1          
          psih  = max( -betah, psih);   %
          uprf  = log(zm2./zm1) - psim;
          uprf  = min(log(zm2./zm1) - psim,2.0*log(zm2./zm1));
          ptprf = min(log(zh2./zh1) - psih,2.0*log(zm2./zm1));
          c_u(lmo>=0.0)  = kv./ uprf;
          c_pt(lmo>=0.0) = kv./ (ptprf * prantl01 );
else %zm1为动力学粗糙度的时候
%#######################################################################
%     Unstable case
%#######################################################################
          rr1=lmo<0.0;
          lmo1=lmo(rr1);
          xx2   = sqrt(sqrt(1-gammam* zm2./lmo1));   
          xx1   = sqrt(sqrt(1-gammam* zm1(rr1)./lmo1));
          psim  = 2 * log((1+xx2)./(1+xx1))+ log((1+xx2.*xx2)./(1+xx1.*xx1))- 2*atan(xx2) + 2 * atan(xx1);
          yy2   = sqrt(1-gammah* zh2./lmo1);  
          yy1   = sqrt(1-gammah* zh1(rr1)./lmo1);
          psih  = 2 * log((1+yy2)./(1+yy1)); 
          uprf  = max(log(zm2./zm1(rr1))- psim,0.50*log(zm2./zm1(rr1)));
          ptprf = max(log(zh2./zh1(rr1)) - psih,0.33*log(zm2./zm1(rr1))); 
          c_u(lmo<0.0)   = kv./ uprf;
          c_pt(lmo<0.0)  = kv./ (ptprf * prantl02 );           
%#######################################################################
%     Stable case:
%#######################################################################
          rr2=lmo>=0.0;
          lmo1=lmo(rr2);
          psim  = -betam * (zm2-zm1(rr2))./lmo1;
          psih  = -betah * (zh2-zh1(rr2))./lmo1;
          psim  = max( -betam, psim);   %实数的对比，fortran用amax1          
          psih  = max( -betah, psih);   %
          uprf  = log(zm2./zm1(rr2)) - psim;
          uprf  = min(log(zm2./zm1(rr2)) - psim,2.0*log(zm2./zm1(rr2)));
          ptprf = min(log(zh2./zh1(rr2)) - psih,2.0*log(zm2./zm1(rr2)));
          c_u(lmo>=0.0)  = kv./ uprf;
          c_pt(lmo>=0.0) = kv./ (ptprf * prantl01 );
end
