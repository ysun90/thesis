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
d.EBB = fitT.EBB'; d.EBBc = d.EBB./sqrt(d.Lepsi.*d.f_plateau);
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
txt={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};

%% WLF vs. nu_eff at different rho in OH 
% for PoP
r=[-.4 0 .4];
dr=0.05;
q=[3 6];
EBB_max = 0.6;
WLF_max = 16;
figure('Pos',[80   200   850   400],'Color','w');  
for jj=1:3
subplot(1,3,jj)
ind=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(1)&d.qpsi<q(2)&...
  d.PLH<0.1 & d.PECRH<0.1 & d.PICRH<0.1 &...
  d.Zeff>1 & d.Zeff<4 &d.EBB<=EBB_max &d.WLF<WLF_max;
hold on
p1=plot(d.nu_eff_l(ind&indLOC),d.WLF(ind&indLOC),'r^','MarkerSize',4);
p2=plot(d.nu_eff_l(ind&indSOC),d.WLF(ind&indSOC),'bs','MarkerSize',4);
p3=plot(d.nu_eff_l(ind&~indLOC&~indSOC),d.WLF(ind&~indLOC&~indSOC),'k.');


% % Plot median
% %N = 200;
% if jj==1
% [xm,ym] = cal_profile_log(d.nu_eff_l(ind),d.WLF(ind),0.3,2,0.1);
% pm1=plot(xm,smooth(ym),'k-.','LineWidth',3,'MarkerSize',16');
% elseif jj==2
% [xm,ym] = cal_profile_log(d.nu_eff_l(ind),d.WLF(ind),0.05,1.2,0.1);
% pm2=plot(xm,smooth(ym),'k-.','LineWidth',3,'MarkerSize',16');
% end

hold off

ax=gca;
set(ax,'FontSize',16,'LineWidth',1.5,'XLim',[0.03 5],'YLim',[2 18],...
   'TickDir','in','XScale','log','TickLength',[.025 .025])
 axis fill
 axis square
box on
%grid on
%ax.Position(1)=0.1;
ax.Position(3)=0.27;
%ax.Position(4)=0.25;
%text(0.035,0.8,txt{ii},'FontSize',18)
%if jj==3
%text(8,0.9,sprintf('%d < q_{\\psi} < %d',q(ii),q(ii+1)),...
 % 'FontSize',16,'Rotation',270);
%end
text(0.035,17,txt{jj},'FontSize',16)
%text(0.035,5,['E_{BB} < ',num2str(EBB_max)],'FontSize',14)
%text(0.035,0.4*3,sprintf('\\rho = %.1f',r(jj)),'FontSize',14);
ax.XTick = [0.05 0.3 1 2 5];  
if jj==1; ax.YTick = [0 5 10 15 20]; else ax.YTickLabel=[]; end
if jj==2; xlabel('\nu_{eff}'); 
  lg=legend([p1,p2,p3],'LOC','SOC','Transition','Location','best');
  lg.FontSize=16;
  lg.Box='on';
  lg.Orientation = 'horizontal';
  lg.Position = [0.3666 0.9049 0.3529 0.0787];
end
if jj==1; ylabel('W_{LF} [kHz]'); end
  if jj==1; title('HFS (\rho = - 0.4)'); end
  if jj==2; title('Center (\rho = 0)'); end
  if jj==3; title('LFS (\rho = 0.4)'); end
end

%% WLF vs. nu_eff at different rho in L-mode
% for PoP
r=[-.4 0 .4];
dr=0.05;
q=[3 6];
P=[1.5 3.5];
EBB_max=0.6;WLF_max=16;
figure('Pos',[80   200   850   400],'Color','w');  
for jj=1:3
subplot(1,3,jj)

ind=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(1)&d.qpsi<q(2)&...
  d.PLH>=P(1) & d.PLH<=P(2) & d.PECRH<0.1 & d.PICRH<0.1 &...
  d.EBB<=EBB_max & d.WLF<=WLF_max;
ind2=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(1)&d.qpsi<q(2)&...
 d.PICRH>=P(1) & d.PICRH<=P(2) & d.PECRH<0.1 & d.PLH<0.1 &...
 d.EBB<=EBB_max  & d.WLF<=WLF_max;

p1=plot(d.nu_eff_l_SL(ind),d.WLF(ind),'r^','MarkerSize',4);hold on
p2=plot(d.nu_eff_l_SL(ind2),d.WLF(ind2),'bs','MarkerSize',4);

% % Plot median
% %N = 200;
% if jj==1
% [xm,ym] = cal_profile_log(d.nu_eff_l_SL(ind),d.WLF(ind),0.1,1,0.1);
% pm1=plot(xm,smooth(ym),'k-.','LineWidth',3,'MarkerSize',16');
% % elseif jj==2
% % [xm,ym] = cal_profile_log(d.nu_eff_l(ind),d.EBB(ind),0.05,1.4,0.1);
% % pm2=plot(xm,smooth(ym),'k-.','LineWidth',3,'MarkerSize',16');
% end

ax=gca;
set(ax,'FontSize',16,'LineWidth',1.5,'XLim',[0.03 5],'YLim',[0 20],...
   'TickDir','in','XScale','log','TickLength',[.025 .025])
 axis fill
 axis square
box on
%grid on
%ax.Position(1)=0.1;
ax.Position(3)=0.27;
%ax.Position(4)=0.25;
%text(0.035,0.8,txt{ii},'FontSize',18)
%if jj==3
%text(8,0.9,sprintf('%d < q_{\\psi} < %d',q(ii),q(ii+1)),...
 % 'FontSize',16,'Rotation',270);
%end
text(0.035,18,txt{jj},'FontSize',16)
text(0.035,5,['E_{BB}<',num2str(EBB_max)],'Fontsize',14)
%text(0.035,0.4*3,sprintf('\\rho = %.1f',r(jj)),'FontSize',14);
ax.XTick = [0.05 0.3 1 2 5];  
if jj==1; ax.YTick = [0 5 10 15 20]; else ax.YTickLabel=[]; end
if jj==2; xlabel('\nu_{eff}'); 
lg=legend([p1,p2],'LH','ICRH','Location','best');
lg.FontSize=16;
lg.Box='on';
lg.Orientation = 'horizontal';
lg.Position = [0.3666 0.9049 0.3529 0.0787];
end
if jj==1; ylabel('W_{LF} [kHz]'); end
  if jj==1; title('HFS (\rho = - 0.4)'); end
  if jj==2; title('Center (\rho = 0)'); end
  if jj==3; title('LFS (\rho = 0.4)'); end
end