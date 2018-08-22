%%
% water = 0; evergreen needleleaf forest = 1;evergreen broadleaf forest = 2;deciduous needleleaf forest = 3;
% deciduous broadleaf forest = 4;mixed forests = 5;closed shrubland = 6; open shrublands = 7;
% woody savannas = 8; savannas = 9; grasslands = 10; permanent wetlands = 11; croplands = 12;
% urban and built-up = 13; cropland/natural vegetation mosaic = 14; snow and ice = 15;
% barren or sparsely vegetated = 16; unclassified = 255;
clc
clear all
constants1
global c1;global c2;global c3;global Ct;global hs1;global Cd;global fh;global A1;global A2;global A3;global xi_m;
global Hi_PBL;global beta0;global zs2hc;global dth;global phiuc; global hcp;global stp;global hs1;global As;
phiuc   =   2;
beta0   =   0.28;
zs2hc   =   3.45;
dth     =   0.7;
fh      =   0;
phiuc   =   3;
Cd      =   0.2;        % Foliage drag coefficient
Ct      =   0.01;       % Heat transfer coefficient
c1      =   0.320;                                           % model constants (Massman 1997)
c2      =   0.264;                                           % model constants (Massman 1997)
c3      =   15.1;                                            % model constants (Massman 1997)
Hi_PBL  =   1000;
hs1     =   0.004; % momentum roughness parameter (0.009 ~ 0.024)(Su et al., 1997, IJRS, p.2105-2124.),default value is 0.0012
kbc     =   5;
kbs     =   5;
A2      =   -5;
As      =   0.5;
cc={'Speud.mat';'North_Carolina_Loblolly_Pine.mat'; 'GLEES.mat';'Marys_River_Fir.mat';'Niwot_Ridge.mat';
    'UMBS.mat';'Chestnut_Ridge_DBF.mat';'Missouri_Ozark.mat';'Morgan_Monroe.mat';'Walker_Branch.mat';'Willow_Creek.mat';
    'Audubon.mat';'Metolius_Young_Pine_Burn.mat';'US_SRM.mat';
    'Gingin.mat';'Calperum.mat';'Ti_Tree_East.mat';'Flagstaff_Managed_Forest.mat';
    'Linzhi.mat';'Brookings.mat';'Cottonwood.mat';'Kendall_Grassland.mat';
    'Riggs_Creek.mat';'Daly_Pasture.mat';'Arcturus.mat';
    'Qolangma.mat';'namco.mat';'maqu.mat'};
for st=6:17
    load(['./',cell2mat(cc(st))]);
    beta =  zeros(size(LAI)); % zeros represent leaf off;, one represent leaf on;
    delh =  0.0001;
    z01  =  0:delh:1;
    nz   =  length(z01);
    stp  =  nz;
    z0m  =  LAI;
    d0   =  LAI;
    z0h  =  LAI;
    for i=1:length(LAI)
        [hcp]     =  leaflength(LCT(i));
        [A2, As, lmdars,Cd]     =  MassmanPa(LCT(i));
        if length(Zref)>1
            Zref  =  Zref(1);
        end
        [z0m(i), d0(i), z0h(i),kB_1(i),kB1_s(i), kB1_c(i),kB1_m(i)]=  kb_1(t1(i),LCT(i),fc(i),NDVI(i),LAI(i),hc(i),Zref,Uref(i),Pref(i),Tref_K(i),LST_K(i),qaref(i),kbc,kbs);
    end
    [Rn,G0,H,LE,ustar] =   EnergyBalance_mapRn(d0, z0m, z0h,LCT,hc, fc,LAI, ...,
        albedo,emissivity,SWd,LWd, LST_K,Hi_PBL, Zref, Tref_K, Uref, Earef,qaref, Pref, Pref,0);
    %---------------------------------------
    H(SWd<10)=NaN;  % valid for daytitme
    %---------------------------------------
    id   =   isnan(tkh) | isnan(H);
    p= polyfit(tkh(~id),H(~id),1);
    RMSEh = sqrt(mean(((tkh(~id)-H(~id)).^2)));
    Re=mean(tkh(~id)-H(~id));
    Cr=corrcoef(tkh(~id),H(~id));
    p1(st)=p(2);p2(st)=p(1);rmse(st)=RMSEh;cr(st)=Cr(2,1);re(1)=Re;
    
    id   =   isnan(tkustar) | isnan(ustar);
    p= polyfit(tkustar(~id),ustar(~id),1);
    RMSEh = sqrt(mean(((tkustar(~id)-ustar(~id)).^2)));
    Re=mean(tkustar(~id)-ustar(~id));
    Cr=corrcoef(tkustar(~id),ustar(~id));
    p11(st)=p(2);p21(st)=p(1);rmse1(st)=RMSEh;cr1(st)=Cr(2,1);re1(1)=Re;
    
    id   =   isnan(tkle) | isnan(LE);
    p= polyfit(tkle(~id),LE(~id),1);
    RMSEh = sqrt(mean(((tkle(~id)-LE(~id)).^2)));
    Re=mean(tkle(~id)-LE(~id));
    Cr=corrcoef(tkle(~id),LE(~id));
    p12(st)=p(2);p22(st)=p(1);rmse2(st)=RMSEh;cr2(st)=Cr(2,1);re2(1)=Re;
end
%%
figure;
h1=subplot(3,1,1);
set(h1,'Unit','Normalized','Position',[0.08,0.74,0.85,0.22]);
% load ????????????????????????_new_ct_massman.mat
plot(p2,'b-','linewidth',1.5);
ylabel('Slope','FontSize',14)
ylim([0 2])
xlim([0.5 28.5])
title('(a)','FontSize',14)
hold on;
plot([0 30],[1 1],'k-')
plot([5+0.5 5+0.5],[0 180],'k-')
plot([11+0.5 11+0.5],[-20 180],'k-')
plot([14+0.5 14+0.5],[-20 180],'k-')
plot([18+0.5 18+0.5],[-20 180],'k-')
plot([22+0.5 22+0.5],[-20 180],'k-')
plot([25+0.5 25+0.5],[-20 180],'k-')
text(0.5,0.1,['ENF'],'Fontname','Times New Roman','FontSize',12);
text(5.5,0.1,['DBF'],'Fontname','Times New Roman','FontSize',12);
text(11.5,0.1,['SRB'],'Fontname','Times New Roman','FontSize',12);
text(14.5,0.1,['SAV'],'Fontname','Times New Roman','FontSize',12);
text(18.5,0.1,['GRS'],'Fontname','Times New Roman','FontSize',12);
text(22.5,0.1,['CRP'],'Fontname','Times New Roman','FontSize',12);
text(25.5,0.1,['BSN'],'Fontname','Times New Roman','FontSize',12);
set(gca,'XTickLabel','')
% legend('H','Ustar','LE')
h1=subplot(3,1,2);
set(h1,'Unit','Normalized','Position',[0.08,0.41,0.85,0.22]); 
plot(rmse,'r-','linewidth',1.5);
hold on;
xlim([0.5 28.5])
ylim([0 200])
ylabel('RMSE','FontSize',14)
title('(b)','FontSize',14)
hold on;
plot([5+0.5 5+0.5],[0 180],'k-')
plot([11+0.5 11+0.5],[-20 180],'k-')
plot([14+0.5 14+0.5],[-20 180],'k-')
plot([18+0.5 18+0.5],[-20 180],'k-')
plot([22+0.5 22+0.5],[-20 180],'k-')
plot([25+0.5 25+0.5],[-20 180],'k-')
text(0.5,4.9,['ENF'],'Fontname','Times New Roman','FontSize',12);
text(5.5,4.9,['DBF'],'Fontname','Times New Roman','FontSize',12);
text(11.5,4.9,['SRB'],'Fontname','Times New Roman','FontSize',12);
text(14.5,4.9,['SAV'],'Fontname','Times New Roman','FontSize',12);
text(18.5,4.9,['GRS'],'Fontname','Times New Roman','FontSize',12);
text(22.5,4.9,['CRP'],'Fontname','Times New Roman','FontSize',12);
text(25.5,4.9,['BSN'],'Fontname','Times New Roman','FontSize',12);
set(gca,'XTickLabel','')
h1=subplot(3,1,3);
set(h1,'Unit','Normalized','Position',[0.08,0.06,0.85,0.23]); 
plot(cr,'r','linewidth',1.5);
hold on;
xlim([0.5 28.5])
ylim([0.3 1])
ylabel('{\itr}','FontSize',14)
title('(c)','FontSize',14)
hold on;
plot([5+0.5 5+0.5],[0 180],'k-')
plot([11+0.5 11+0.5],[-20 180],'k-')
plot([14+0.5 14+0.5],[-20 180],'k-')
plot([18+0.5 18+0.5],[-20 180],'k-')
plot([22+0.5 22+0.5],[-20 180],'k-')
plot([25+0.5 25+0.5],[-20 180],'k-')
text(0.5,0.35,['ENF'],'Fontname','Times New Roman','FontSize',12);
text(5.5,0.35,['DBF'],'Fontname','Times New Roman','FontSize',12);
text(11.5,0.35,['SRB'],'Fontname','Times New Roman','FontSize',12);
text(14.5,0.35,['SAV'],'Fontname','Times New Roman','FontSize',12);
text(18.5,0.35,['GRS'],'Fontname','Times New Roman','FontSize',12);
text(22.5,0.35,['CRP'],'Fontname','Times New Roman','FontSize',12);
text(25.5,0.35,['BSN'],'Fontname','Times New Roman','FontSize',12);
set(gca,'XTickLabel','')
