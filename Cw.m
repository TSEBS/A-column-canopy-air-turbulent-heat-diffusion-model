function [cw] = Cw(hpbl, L, z0m, z0h, d0) 
% Bulk Stability function for heat tranfer, eq.(23), (27)
% hpbl: Height of ABL or PBL
% L: The Obukhov length
% z0: Surface roughnes height for momentum
% z0h: Surface roughnes height for height transfer
% 
% The ASL height
% hst = alfa*hpbl, with alfa=0.10~0.15 over moderately rough surfaces, or hst=beta*z0, with beta= 100~150.
alfa                                =   0.12;                              % These are mid values as given by Brutsaert, 1999
beta                                =   125 ;                              % These are mid values as given by Brutsaert, 1999
zhl                                 =   (-z0h ./ L);                       % 分开写可以加快运算速度
I_s                                 =   zhl < 0;                           % STABLE conditions (Brutsaert, 1982, Eq. 4.93,p.84)
cw_s                                =   -7.6 * log(1 -zhl);      
I_u                                 =   zhl >=0;                           % UNSTABLE conditions (Brutsaert 2008, p53)       
abh                                 =   (alfa ./ beta) * hpbl;
I_umr                               =   z0m < abh;                         % for moderately rough terrain (Brutsaert 2008, p53)       
I_uvr                               =   z0m >=abh;                         % for very rough terrain (Brutsaert 2008, p53)
tp1                                 =   alfa .* (hpbl-d0)./L;
tp2                                 =   z0h./L;
psih1                               =   PSIh(tp1);
psih2                               =   PSIh(tp2);  
cw_umr                              =   psih1 - psih2 - log(alfa);
tp1                                 =   beta * (z0m  )./L;
tp2                                 =   z0h./L;
tp3                                 =   (hpbl-d0)./(beta*z0m);
psih1                               =   PSIh(tp1);
psih2                               =   PSIh(tp2);                           
cw_uvr                              =   psih1 - psih2 + log(tp3);
tp                                  =   I_umr.*cw_umr + I_uvr.*cw_uvr;
% cw_u                                =   max(tp ,0);
cw_u                                =   tp;  % 代替ｍａｘ的计算，运算速度
cw_u(cw_u<0)                        =   0;
cw                                  =   I_s.*cw_s + I_u.*cw_u;
return