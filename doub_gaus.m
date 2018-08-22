function  [xi_m,sd_u,sd_l]   =  doub_gaus(LCT,beta,varargin)
global ksim1;global sdu1;global sdl1;
% xi_m, sd_u, sd_l
% Massman 2016
% Summer
%              water   ENF    EBF   DNF    DBF    MF    SRB    SRB    SAV    SAV    GRS    WET    CRP    URB     CRP   snow  bare   
ksim1    =   [0.200; 0.600; 0.600; 0.550; 0.550; 0.600; 0.95; 0.95; 0.400; 0.400; 0.990; 0.200; 0.720; 0.600;  0.720; 0.600; 0.900; 0.700];
sdu1     =   [0.010; 0.180; 0.180; 0.400; 0.400; 0.150; 0.35; 0.35; 0.150; 0.150; 0.550; 0.001; 0.010; 0.050;  0.010; 0.001; 0.140; 0.040];
sdl1     =   [0.001; 0.060; 0.060; 0.300; 0.300; 0.060; 0.001; 0.95; 0.050; 0.050; 0.030; 0.001; 0.001; 0.050;  0.0010; 0.001;0.001; 0.001];
% Winter
%               water  ENF    EBF   DNF    DBF     MF    SRB    SRB    SAV    SAV    GRS    WET    CRP    URB     CRP    snow  bare   
ksim2    =   [0.200; 0.600; 0.600; 0.600; 0.600; 0.700; 0.95; 0.95; 0.450; 0.450; 0.990; 0.200; 0.720; 0.600;  0.720; 0.600; 0.700; 0.700];
sdu2     =   [0.010; 0.080; 0.080; 0.080; 0.080; 0.150; 0.95; 0.95; 0.550; 0.550; 0.550; 0.001; 0.010; 0.050;  0.010; 0.001; 0.040; 0.040];
sdl2     =   [0.001; 0.060; 0.060; 0.060; 0.060; 0.060; 0.95; 0.95; 0.100; 0.050; 0.551; 0.001; 0.001; 0.050;  0.001; 0.001; 0.001; 0.001];
ksim1    =   single(ksim1); ksim2   =   single(ksim2);
sdu1     =   single(sdu1);  sdu2    =   single(sdu2);
sdl1     =   single(sdl1);  sdl2    =   single(sdl2);
sd_u     =   LCT;
sd_l     =   LCT;
xi_m     =   LCT;
if nargin == 1
    for i=0:17
        sd_u(LCT==i)  = sdu1(i+1);  % 
        sd_l(LCT==i)  = sdl1(i+1);  %
        xi_m(LCT==i)  = ksim1(i+1); %
    end
elseif nargin == 2
    for i=0:17
        sd_u(LCT==i&beta==1)  = sdu1(i+1); % 
        sd_l(LCT==i&beta==1)  = sdl1(i+1);%
        xi_m(LCT==i&beta==1)  = ksim1(i+1); %
        sd_u(LCT==i&beta==0)  = sdu2(i+1); % 
        sd_l(LCT==i&beta==0)  = sdl2(i+1);%
        xi_m(LCT==i&beta==0)  = ksim2(i+1); %
    end
end