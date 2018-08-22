function [lmo]=MOlength(zm,zh,z0m,z0h,wspd,ptsfc,pt1)
%    MOlength(zm,zh,z0m,z0h,wspd,ptsfc, pt1)
%    Compute lmo for both unstable and stable surface layer 
%    based on Monin-Obukhov similarity theory. The univeral profile form 
%    was proposed Hogstrom (1996) and the anlytical solution was 
%    developed by Yang (2001)
%    REFERENCE: Similarity analytical solution:Yang, K., Tamai, N. and Koike, T., 2001: Analytical Solution 
%         of Surface Layer Similarity Equations, J. Appl. Meteorol. 40, 
%         1647-1653. Profile form: Hogstrom (1996,Boundary-Layer Meteorol.)
g         =  9.8; 
prantl01  =  1.0;  % Turbulent Prandtl number for stable case
prantl02  =  0.95; % Turbulent Prandtl number for unstable case
gammam    =  19.0;
gammah    =  11.6;
betam     =  5.3;
betah     =  8.0;
%
[m,n]     =  size(wspd);
lmo(1:m,1:n)=NaN;
%
bulkri    =  (g./pt1).*(pt1-ptsfc).*(zm-z0m)./(wspd.^2);
rr1       =  bulkri<0.0;
rr2       =  bulkri>=0.0;
%#######################################################################
%    Unstable case: calculated by an anproximate analytical solution
%    proposed by Yang (2001,JAM)%
%#######################################################################%
bulkri1   =  max(bulkri(rr1),-10.0);     
d         =  bulkri1./prantl02;
numerator =  d.*(log(zm./z0m(rr1)).^2./log(zh./z0h(rr1))).*(1./(zm-z0m(rr1)));  %转化过程中要主意平方和其他算符的顺序 
a         =  log(-d);
b         =  log(log(zm./z0m(rr1)));
c         =  log(log(zh./z0h(rr1)));		 
p         =  0.03728-0.093143*a-0.24069*b+0.30616*c+0.017131*a.*a+0.037666*a.*b -0.084598*b.*b-0.016498*a.*c+0.1828*b.*c-0.12587*c.*c ;
p         =  max(0.0,p) ;     
coef      =  d*gammam^2/8/gammah.*(zm-z0m(rr1))./(zh-z0h(rr1));  
lmo(rr1)  =  numerator./(1-coef.*p);
%#######################################################################
%    Stable case: calculated by the exact analytical solution
%    proposed by Yang (2001,JAM)%
%#######################################################################
bulkri1  =  min(bulkri(rr2),min(0.2,prantl01*betah*(1-z0h(rr2)/zh)/betam^2./(1-z0m(rr2)/zm)-0.05));  
d        =  bulkri1/prantl01;
a        =  d*betam^2.*(zm-z0m(rr2))-betah*(zh-z0h(rr2));
b        =  2*d*betam.*log(zm./z0m(rr2))-log(zh./z0h(rr2));
c        =  d.*log(zm./z0m(rr2)).^2./(zm-z0m(rr2));
lmo(rr2) =  (-b-sqrt(b.*b-4*a.*c))./(2*a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmo(lmo>0&lmo<1.0e-6)=1.0e-6; %设定lmo的范围，大于零小于0.000006时，设定为0.000006
lmo(lmo<0&lmo>-1.0e-6)=-1.0e-6; %设定lmo的范围，小于零大于-0.000006时，设定为-0.000006
lmo=1./lmo;
