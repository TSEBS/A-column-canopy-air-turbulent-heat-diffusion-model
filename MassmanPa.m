function  [A2,as,lmdars,cd,varargout]   =  MassmanPa(LCT,varargin)
% Reference: Chen et al. 2017 AFM submitted.
%        water  ENF  EBF  DNF  DBF   MF  SRB   SRB   SAV   SAV  GRS  WET  CRP   URB  CRP   snow  bare   
a2    =   [-5.; -5.; -5.; -5.; -5.; -5.; -5.;  -5.;  -5.;  -5.; -5.; -5.; -5.;  -5.;  -5.;  -5.;  -5.; -5. ];
AS    =   [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.7;  0.7;  0.5;  0.5; 0.7; 0.5; 0.5;  0.5;  0.5;  0.5;  0.7; 0.5 ];
% AS    =   [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;  0.0;  0.0;  0.0; 0.0; 0.0;  0.0; 0.0;  0.0;  0.0;  0.0; 0.0];
% AS    =   [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1;  0.1;  0.1;  0.1; 0.1; 0.1;  0.1; 0.1;  0.1;  0.1;  0.1; 0.1];
ls    =   [  1; 15.; 15.; 15.; 15.; 15.; 11.;  1.1;  1.1;  1.1; 1.1; 1.1;  1.1; 1.1;  1.1;  1.1;  1.1; 1.1];
Cd    =   [0.2; 0.2; 0.2; 0.2; 0.2; 0.2;0.11; 0.15;  0.2;  0.2; 0.2; 0.2;  0.2; 0.2;  0.2;  0.2;  0.15; 0.2];

A2  = LCT;
as  = LCT;
lmdars= LCT;
cd=LCT;
for i=0:17
    as(LCT==i)  = AS(i+1); % 
    A2(LCT==i)  = a2(i+1); % 
    lmdars(LCT==i)  = ls(i+1); % 
    cd(LCT==i)  = Cd(i+1); % 
end

return