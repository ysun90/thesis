%% Turbulence database by Taylor model
clear
% Reduced database
load('DB.mat'); d=DB;
db = load('Database_v10_plasma_temp.mat');
load('db_Zeff.mat')
load('db_Zeff_TAN.mat')
load('db_Tmoy.mat')
d.Tmoy=Tmoy;
d.rhoc = db.reflec.rhoc;
d.Nec=db.plasma.Nec;
d.Tec=db.plasma.Tec;
% Averaged linear density
d.Nla = d.Nl./(2*d.a); % [10^19 m-3]
% Normalized radius
load('rhoc.mat'); d.rhoc=r;
load('db_Nelc.mat')
load('db_icrh.mat');
d.PICRH=db_pos_icrh.p_icrh_total;
% Taylor model
fitT = load('T358874.mat');
fitG = load('GG358874.mat');
d.EBB = fitT.EBB'; 
d.WBB=fitT.SDBB';
d.WLF=fitT.sigmaLF';
d.betaBB=fitG.betaBB';
d.k2DBB=fitT.k2DBB'; %(k^2*u^2*tau)
d.tauBB=fitT.tauBB';d.kutauBB=sqrt(d.k2DBB.*d.tauBB);
d.ELF = fitT.ELF';d.sigmaLF=fitT.sigmaLF';d.muLF = fitT.muLF';
d.ECS=fitT.ECS';d.sigmaCS=fitT.sigmaCS';
d.ES0=fitT.ES0';
d.muCS = fitT.muCS';
d.muBB = fitT.muBB';
d.AmpBB=fitT.AmpBB';
d.Noise=fitT.Noise';
d.FVAL=fitT.FVAL';
d.Nelc=Nelc_robust;
d.Zeff=Zeff_TAN;
d.Zeff_SL=1+7*d.Prad./((d.Nl*1e19*1e-20).^2.*d.a.*pi^2.*d.Rp*4);
%d.Zeff=d.Zeff_pred;
d.nu_eff_l_Z2=0.1*2*d.Rp.*(d.Nec/1e19)./(d.Tec.^2);
d.nu_eff_l=0.1*d.Zeff.*d.Rp.*(d.Nec/1e19)./(d.Tec.^2);
d.nu_eff_l_SL=0.1*d.Zeff_SL.*d.Rp.*(d.Nec/1e19)./(d.Tec.^2);
%d.nu_eff_l=0.1*d.Zeff.*d.Rp.*(d.Nmoy/1e19)./(d.Tmoy.^2);
d.nu_ei_Z2=0.1*2.*d.Rp.*(d.Nec/1e19)./(d.Tec.^1.5);
d.nu_ei=0.1*d.Zeff.*d.Rp.*(d.Nec/1e19)./(d.Tec.^1.5);
clear DB reflec fitT r Nelc Nelc_robust db 
%d.EBB=d.EBB.*d.ES0;

% Filter data
 %indOK = d.FVAL<0.002 & lg(d.EBB./(d.Noise*1025))>0 & d.EquiTag==1 & ...
    % abs(d.muBB)<50;
  %indOK = d.FVAL<0.002 & lg(d.EBB./(d.Noise*1025))>0 & d.EquiTag==1 & ...
   % abs(d.muBB)<50 & d.Tec>0;
 indOK = d.FVAL<0.002 & lg(d.AmpBB./d.Noise)>25 & d.EquiTag==1 & ...
    abs(d.muBB)<50 & d.Nla>0 & d.f_plateau<1e3;
 % indOK = d.FVAL<0.002 & lg(d.AmpBB./d.Noise)>25 & d.EquiTag==1 & ...
  %  abs(d.muBB)<50 & d.Nla>0 & d.f_plateau<1e3 & (d.Zeff>=1 & d.Zeff<=4);
% Ohmic 
indOH = (d.PICRH + d.PECRH + d.PLH)<0.1;
indLOC = d.Nla<2.6*d.Ip*0.9 & indOH; 
indSOC = d.Nla>2.6*d.Ip*1.1 & indOH;
indLS = {indLOC,indSOC};
% PICH
indICRH1 = d.PICRH>0.5 & d.PICRH<1.5 & d.PECRH<0.1 & d.PLH<0.1;
indICRH2 = d.PICRH>1.5 & d.PICRH<2.5 & d.PECRH<0.1 & d.PLH<0.1;
indICRH3 = d.PICRH>2.5 & d.PICRH<3.5 & d.PECRH<0.1 & d.PLH<0.1;
indICRH4 = d.PICRH>3.5 & d.PICRH<inf & d.PECRH<0.1 & d.PLH<0.1;
indICRH5 = d.PICRH>0.5 & d.PICRH<inf & d.PECRH<0.1 & d.PLH<0.1;
indICRH = {indICRH1,indICRH2,indICRH3,indICRH4,indICRH5};
% PLH
indLH1 = d.PLH>0.5 & d.PLH<1.5 & d.PECRH<0.1 & d.PICRH<0.1;
indLH2 = d.PLH>1.5 & d.PLH<2.5 & d.PECRH<0.1 & d.PICRH<0.1;
indLH3 = d.PLH>2.5 & d.PLH<3.5 & d.PECRH<0.1 & d.PICRH<0.1;
indLH4 = d.PLH>3.5 & d.PLH<inf & d.PECRH<0.1 & d.PICRH<0.1;
indLH5 = d.PLH>0.5 & d.PLH<inf & d.PECRH<0.1 & d.PICRH<0.1;
indLH = {indLH1,indLH2,indLH3,indLH4,indLH5};
% PECRH
indECRH=d.PLH<0.5 & d.PECRH>0.1 & d.PICRH<0.5;
indECRH = {indECRH,indECRH,indECRH,indECRH};
% Safety factor 
q = [3 4 5 6 10];
indq = cell(1,length(q)-1);
for i = 1:length(indq)
    indq{i} = d.qpsi>q(i) & d.qpsi<q(i+1);
end
% figure label
txt={'(a)','(b)','(c)','(d)'};

%% Radial profiles of WBB and betaBB
r=-.95:.05:.55; dr=.05;
q=[3 4 5 6];
P=[1.5 3.5];
plotLH={'r-','r--','r:'};plotICRH={'b-','b--','b:'};
figure('Pos',[96         137        1187         503],'Color','w');

subplot(211)
for ii=1:3
for jj=1:length(r)
ind=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
  d.PLH>=P(1) & d.PLH<=P(2) & d.PECRH<0.1 & d.PICRH<0.1;
WBB_mean_LH(jj)=median(d.WBB(ind));
ind2=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
 d.PICRH>=P(1) & d.PICRH<=P(2) & d.PECRH<0.1 & d.PLH<0.1;
  if length(d.WBB(ind2))<5
   WBB_mean_ICRH(jj)=NaN; 
  else
   WBB_mean_ICRH(jj)=median(d.WBB(ind2));  
  end
end
%plot(r,WBB_mean_LH,plotLH{ii},'LineWidth',2);hold on
if ii==1; plot(r(r<.3),WBB_mean_ICRH(r<.3),plotICRH{ii},'LineWidth',2); hold on
else
plot(r,WBB_mean_ICRH,plotICRH{ii},'LineWidth',2); hold on
end
end
xlabel('\rho');ylabel('W_{BB}');ylim([50 200])
legend('3 < q < 4', '4 < q < 5','5 < q < 6');
% legend('LH, 3 < q < 4', 'ICRH, 3 < q < 4','LH, 4 < q < 5',...
%   'ICRH, 4 < q < 5', 'LH, 5 < q < 6','ICRH, 5 < q < 6');

subplot(212)
for ii=1:3
  
for jj=1:length(r)
ind=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
  d.PLH>=P(1) & d.PLH<=P(2) & d.PECRH<0.1 & d.PICRH<0.1;
betaBB_mean_LH(jj)=median(d.betaBB(ind));
ind2=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
 d.PICRH>=P(1) & d.PICRH<=P(2) & d.PECRH<0.1 & d.PLH<0.1;
  if length(d.betaBB(ind2))<5
   betaBB_mean_ICRH(jj)=NaN; 
  else
   betaBB_mean_ICRH(jj)=median(d.betaBB(ind2));  
  end
end
%plot(r,betaBB_mean_LH,plotLH{ii},'LineWidth',2);hold on;
if ii==1; plot(r(r<.3),betaBB_mean_ICRH(r<.3),plotICRH{ii},'LineWidth',2); hold on
else
plot(r,betaBB_mean_ICRH,plotICRH{ii},'LineWidth',2); hold on
end
end
xlabel('\rho');ylabel('\beta_{BB}');ylim([1 2])
legend('3 < q < 4', '4 < q < 5','5 < q < 6');