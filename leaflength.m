function [LL]  =    leaflength(LCT)
% water = 0; evergreen needleleaf forest = 1;evergreen broadleaf forest = 2;deciduous needleleaf forest = 3;
% deciduous broadleaf forest = 4;mixed forests = 5;closed shrubland = 6; open shrublands = 7;
% woody savannas = 8; savannas = 9; grasslands = 10; permanent wetlands = 11; croplands = 12;
% urban and built-up = 13; cropland/natural vegetation mosaic = 14; snow and ice = 15;
% barren or sparsely vegetated = 16; unclassified = 255;
%       water  ENF    EBF   DNF     DBF     MF    SRB   SRB    SAV    SAV    GRS    WET     CRP    URB    CRP    snow   bare   
leafL = [1.000; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000; 1.000;  1.000; 1.000; 1.000; 1.000];
LL  = LCT;
for i=0:17
    LL(LCT==i)  = leafL(i+1)/100; % 
end
% LL  =  0.005;
