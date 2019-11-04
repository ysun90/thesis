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

%% EBB/WBB/betaBB vs. nu_eff at different q and rho in OH
% for PoP
r=[-.4 0 .4];
dr=0.05;
q=[3 4 5 6];
%Nl=[2 8]*1e19;
figure('Pos',[183   121   912   503],'Color','w');  
for ii=1:3
for jj=1:3
subplot(3,3,sub2ind([3,3],jj,ii))
ind=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
  d.PLH<0.1 & d.PECRH<0.1 & d.PICRH<0.1 & d.Zeff>1 & d.Zeff<4;
hold on
plot(d.nu_eff_l(ind&indLOC),d.betaBB(ind&indLOC),'r^','MarkerSize',4);
plot(d.nu_eff_l(ind&indSOC),d.betaBB(ind&indSOC),'bo','MarkerSize',4);
plot(d.nu_eff_l(ind&~indLOC&~indSOC),d.betaBB(ind&~indLOC&~indSOC),'k.');
hold off

ax=gca;
set(ax,'FontSize',12,'LineWidth',1.5,'XLim',[0.03 5],'YLim',[0.2 3.1],...
   'TickDir','out','XScale','log')
%box on
%grid on
%ax.Position(1)=0.15;
ax.Position(3)=0.23;
ax.Position(4)=0.25;
%lg=legend('LOC','Transition','SOC','Location','southeast');
%lg.FontSize=14;
%lg.Orientation = 'horizontal';
%text(0.035,0.8,txt{ii},'FontSize',18)
text(0.035,0.8*3,sprintf('%d < q_{\\psi} < %d',q(ii),q(ii+1)),'FontSize',14);
text(0.035,0.4*3,sprintf('\\rho = %.1f',r(jj)),'FontSize',14);
if ii==3; ax.XTick = [0.05 0.3 1 2 5]; else ax.XTickLabel = []; end
if ii==3 & jj==2; xlabel('\nu_{eff}'); end
if ii==2 & jj==1; ylabel('\beta_{BB}'); end
end
end

%% EBB/WBB/betaBB vs. nu_eff at different q and rho in L-mode
% for PoP
r=[-.4 0 .4];
dr=0.05;
q=[3 4 5 6];
P=[1.5 3.5];
figure('Pos',[183   121   912   503],'Color','w');  
for ii=1:3
for jj=1:3
subplot(3,3,sub2ind([3,3],jj,ii))

ind=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
  d.PLH>=P(1) & d.PLH<=P(2) & d.PECRH<0.1 & d.PICRH<0.1;
ind2=indOK&d.rhoc>r(jj)-dr&d.rhoc<r(jj)+dr&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
 d.PICRH>=P(1) & d.PICRH<=P(2) & d.PECRH<0.1 & d.PLH<0.1;

plot(d.nu_eff_l_SL(ind),d.WBB(ind),'r^','MarkerSize',4);hold on
plot(d.nu_eff_l_SL(ind2),d.WBB(ind2),'bo','MarkerSize',4);

ax=gca;
set(ax,'FontSize',12,'LineWidth',1.5,'XLim',[0.03 5],'YLim',[20 180],...
   'XScale','log','XTickLabelMode','auto','XTick',[0.05 0.3 1 2],...
   'TickLength',[.02 .025])
%box on
grid on
%ax.Position(1)=0.15;
ax.Position(3)=0.23;
ax.Position(4)=0.25;
%lg=legend('LOC','Transition','SOC','Location','southeast');
%lg.FontSize=14;
%lg.Orientation = 'horizontal';
%text(0.035,0.8,txt{ii},'FontSize',18)
text(0.035,0.8*150,sprintf('%d < q_{\\psi} < %d',q(ii),q(ii+1)),'FontSize',14);
text(0.035,0.4*150,sprintf('\\rho = %.1f',r(jj)),'FontSize',14);
if ii==3 & jj==2; xlabel('\nu_{eff}'); end
if ii==2 & jj==1; ylabel('W_{BB}'); end
end
end


%% E_BB vs. nu_eff in different r, Nl, with OH at fixed q
% for PoP
r=[-0.55 -0.45 -0.15 -0.05 0.45 0.55];
q=[4 5];
Nl=[2 8]*1e19;
%P=[1.5 3.5];
%P=[0.5 1.5 2.5 3.5 10];

for kk=1:length(Nl)-1
figure('Pos',[60    85   558   590],'Color','w');  
for ii=1:length(r)/2
subplot(length(r)/2,1,ii)
ind=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
  d.PLH<0.1 & d.PECRH<0.1 & d.PICRH<0.1 &...
   d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
hold on
plot(d.nu_eff_l(ind&indLOC),d.EBB(ind&indLOC),'r^');
plot(d.nu_eff_l(ind & ~indLOC & ~indSOC),d.EBB(ind & ~indLOC & ~indSOC),'kd');
plot(d.nu_eff_l(ind&indSOC),d.EBB(ind&indSOC),'bo');
hold off

ax=gca;
set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[0.05 10],'YLim',[0 1],...
   'TickDir','out','XScale','log')
box on
%grid on
ax.Position(1)=0.15;
ax.Position(4)=0.25;
ylabel('E_{BB}') 
lg=legend('LOC','Transition','SOC','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(0.035,0.8,txt{ii},'FontSize',18)
text(0.035,0.4,sprintf('%.2f < r < %.2f',r(ii*2-1),r(ii*2)),'FontSize',14);
text(0.035,0.2,'4 < q_{\psi} < 5','FontSize',14);
if ii==3; ax.XTick = [0.05 0.1 0.2 0.5 1 2 5 10]; xlabel('\nu_{eff}');
else ax.XTickLabel = [];
end
end
end
% for i=1:48; figure(i); saveas(gcf,['f',num2str(i)],'png'); end

%% W_BB vs. nu_eff in different r, Nl, with OH at fixed q
% for PoP
r=[-0.55 -0.45 -0.15 -0.05 0.45 0.55];
q=[4 5];
Nl=[2 8]*1e19;
%P=[1.5 3.5];
%P=[0.5 1.5 2.5 3.5 10];

for kk=1:length(Nl)-1
figure('Pos',[60    85   558   590],'Color','w');  
for ii=1:length(r)/2
subplot(length(r)/2,1,ii)
ind=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
  d.PLH<0.1 & d.PECRH<0.1 & d.PICRH<0.1 &...
   d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
hold on
plot(d.nu_eff_l(ind&indLOC),d.WBB(ind&indLOC),'r^');
plot(d.nu_eff_l(ind & ~indLOC & ~indSOC),d.WBB(ind & ~indLOC & ~indSOC),'kd');
plot(d.nu_eff_l(ind&indSOC),d.WBB(ind&indSOC),'bo');
hold off

ax=gca;
set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[0.05 10],'YLim',[20 180],...
   'TickDir','out','XScale','log')
box on
%grid on
ax.Position(1)=0.15;
ax.Position(4)=0.25;
ylabel('W_{BB}') 
lg=legend('LOC','Transition','SOC','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(0.1,150,txt{ii},'FontSize',18)
text(0.1,100,sprintf('%.2f < r < %.2f',r(ii*2-1),r(ii*2)),'FontSize',14);
text(0.1,50,'4 < q_{\psi} < 5','FontSize',14);
if ii==3; ax.XTick = [0.05 0.1 0.2 0.5 1 2 5 10]; xlabel('\nu_{eff}');
else ax.XTickLabel = [];
end
end
end
% for i=1:48; figure(i); saveas(gcf,['f',num2str(i)],'png'); end

%% beta_BB vs. nu_eff in different r, Nl, with OH at fixed q
% for PoP
r=[-0.55 -0.45 -0.15 -0.05 0.45 0.55];
q=[4 5];
Nl=[2 8]*1e19;
%P=[1.5 3.5];
%P=[0.5 1.5 2.5 3.5 10];

for kk=1:length(Nl)-1
figure('Pos',[60    85   558   590],'Color','w');  
for ii=1:length(r)/2
subplot(length(r)/2,1,ii)
ind=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
  d.PLH<0.1 & d.PECRH<0.1 & d.PICRH<0.1 &...
   d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
hold on
plot(d.nu_eff_l(ind&indLOC),d.betaBB(ind&indLOC),'r^');
plot(d.nu_eff_l(ind & ~indLOC & ~indSOC),d.betaBB(ind & ~indLOC & ~indSOC),'kd');
plot(d.nu_eff_l(ind&indSOC),d.betaBB(ind&indSOC),'bo');
hold off

ax=gca;
set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[0.05 10],'YLim',[0.5 2.5],...
   'TickDir','out','XScale','log')
box on
%grid on
ax.Position(1)=0.15;
ax.Position(4)=0.25;
ylabel('\beta_{BB}') 
lg=legend('LOC','Transition','SOC','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(0.1,150,txt{ii},'FontSize',18)
text(0.1,100,sprintf('%.2f < r < %.2f',r(ii*2-1),r(ii*2)),'FontSize',14);
text(0.1,50,'4 < q_{\psi} < 5','FontSize',14);
if ii==3; ax.XTick = [0.05 0.1 0.2 0.5 1 2 5 10]; xlabel('\nu_{eff}');
else ax.XTickLabel = [];
end
end
end
% for i=1:48; figure(i); saveas(gcf,['f',num2str(i)],'png'); end

%% E_BB vs. nu_eff in different q, Nl, Power with P_ICRH and P_LH
% for PoP
%r=[-0.75 -0.55 -0.35 -0.15 0.05 0.20 0.40 0.60];%r=[-0.15 0.05];
r=[-0.15 0.05];
q=[3 4 5 6];
Nl=[2 8]*1e19;
P=[1.5 3.5];
%P=[0.5 1.5 2.5 3.5 10];

for jj=1:length(P)-1
for kk=1:length(Nl)-1
for ll=1:length(r)-1
figure('Pos',[60    85   558   590],'Color','w');  
for ii=1:length(q)-1
 subplot(length(q)-1,1,ii)
%subplot(121) % LH
ind=indOK&d.rhoc>r(ll)&d.rhoc<r(ll+1)&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
  d.PLH>P(jj) & d.PLH<P(jj+1) & d.PECRH<0.1 & d.PICRH<0.1 &...
   d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
IND{ii}=ind;
plot(d.nu_eff_l(ind),d.EBB(ind),'r^');hold on

%ax=gca;
% set(ax,'FontSize',14,'LineWidth',1.5,'XLim',[0.03 2],'YLim',[0 1],'XScale','log')
%title(sprintf('%d < q < %d',q(ii),q(ii+1)))
% title(sprintf('%d < q < %d, %.2g < Nl < %.2g, %.2g < P < %.2g, %.2g < r/a < %.2g',...
%   q(ii),q(ii+1),Nl(kk),Nl(kk+1),P(jj),P(jj+1),r(1),r(2)))

ind2=indOK&d.rhoc>r(ll)&d.rhoc<r(ll+1)&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
 d.PICRH>P(jj) & d.PICRH<P(jj+1) & d.PECRH<0.1 & d.PLH<0.1 &...
 d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
 plot(d.nu_eff_l(ind2),d.EBB(ind2),'bo');

ax=gca;
 set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[0.03 2],'YLim',[0 1],...
   'TickDir','out','XScale','log')
if ii==3; xlabel('\nu_{eff}'); else ax.XTickLabel = []; end
ax.Position(1)=0.15;
ax.Position(4)=0.25;
ylabel('E_{BB}') 
lg=legend('LH','ICRH','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(0.035,0.8,txt{ii},'FontSize',18)
text(0.035,0.4,sprintf('%d < q_{\\psi} < %d',q(ii),q(ii+1)),'FontSize',14);
text(0.035,0.2,'\rho = - 0.1 \pm 0.05','FontSize',14);
%text(0.035,0.2,sprintf('%.2f < \\rho < %.2f',r(ll),r(ll+1)),'FontSize',14);
if ii==3; ax.XTick = [0.05 0.1 0.2 0.5 1 2]; xlabel('\nu_{eff}');
else ax.XTickLabel = [];end
end
end
end
end

%% E_BB vs. nu_eff in different r, Nl, Power with P_ICRH and P_LH
% for PoP
%r=[-0.75 -0.55 -0.35 -0.15 0.05 0.20 0.40 0.60];%r=[-0.15 0.05];
r=[-0.55 -0.45 -0.15 -0.05 0.45 0.55];
q=[4 5];
Nl=[2 8]*1e19;
P=[1.5 3.5];
%P=[0.5 1.5 2.5 3.5 10];

for jj=1:length(P)-1
for kk=1:length(Nl)-1

figure('Pos',[60    85   558   590],'Color','w');  
for ii=1:length(r)/2
 subplot(length(r)/2,1,ii)
%subplot(121) % LH
ind=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
  d.PLH>P(jj) & d.PLH<P(jj+1) & d.PECRH<0.1 & d.PICRH<0.1 &...
   d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
IND{ii}=ind;
plot(d.nu_eff_l(ind),d.EBB(ind),'r^');hold on

%ax=gca;
% set(ax,'FontSize',14,'LineWidth',1.5,'XLim',[0.03 2],'YLim',[0 1],'XScale','log')
%title(sprintf('%d < q < %d',q(ii),q(ii+1)))
% title(sprintf('%d < q < %d, %.2g < Nl < %.2g, %.2g < P < %.2g, %.2g < r/a < %.2g',...
%   q(ii),q(ii+1),Nl(kk),Nl(kk+1),P(jj),P(jj+1),r(1),r(2)))

ind2=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
 d.PICRH>P(jj) & d.PICRH<P(jj+1) & d.PECRH<0.1 & d.PLH<0.1 &...
 d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
 plot(d.nu_eff_l(ind2),d.EBB(ind2),'bo');

ax=gca;
 set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[0.05 10],'YLim',[0 1],...
   'TickDir','out','XScale','log')
if ii==3; xlabel('\nu_{eff}'); else ax.XTickLabel = []; end
ax.Position(1)=0.15;
ax.Position(4)=0.25;
ylabel('E_{BB}') 
lg=legend('LH','ICRH','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(0.05,0.8,txt{ii},'FontSize',18)
text(0.05,0.4,sprintf('%.2f < r < %.2f',r(ii*2-1),r(ii*2)),'FontSize',14);
text(0.05,0.2,'4 < q_{\psi} < 5','FontSize',14);
%text(0.035,0.2,sprintf('%.2f < \\rho < %.2f',r(ll),r(ll+1)),'FontSize',14);
if ii==3; ax.XTick = [0.05 0.1 0.2 0.5 1 2 5 10]; xlabel('\nu_{eff}');
else ax.XTickLabel = [];end

end
end
end

%% W_BB vs. nu_eff in different r, Nl, Power with P_ICRH and P_LH
% for PoP
%r=[-0.75 -0.55 -0.35 -0.15 0.05 0.20 0.40 0.60];%r=[-0.15 0.05];
r=[-0.55 -0.45 -0.15 -0.05 0.45 0.55];
q=[4 5];
Nl=[2 8]*1e19;
P=[1.5 3.5];
%P=[0.5 1.5 2.5 3.5 10];

for jj=1:length(P)-1
for kk=1:length(Nl)-1

figure('Pos',[60    85   558   590],'Color','w');  
for ii=1:length(r)/2
 subplot(length(r)/2,1,ii)
%subplot(121) % LH
ind=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
  d.PLH>P(jj) & d.PLH<P(jj+1) & d.PECRH<0.1 & d.PICRH<0.1 &...
   d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
plot(d.nu_eff_l(ind),d.WBB(ind),'r^');hold on

ind2=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
 d.PICRH>P(jj) & d.PICRH<P(jj+1) & d.PECRH<0.1 & d.PLH<0.1 &...
 d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
 plot(d.nu_eff_l(ind2),d.WBB(ind2),'bo');

ax=gca;
 set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[0.05 10],'YLim',[20 180],...
   'TickDir','out','XScale','log')
if ii==3; xlabel('\nu_{eff}'); else ax.XTickLabel = []; end
ax.Position(1)=0.15;
ax.Position(4)=0.25;
ylabel('W_{BB}') 
lg=legend('LH','ICRH','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(0.05,0.8,txt{ii},'FontSize',18)
text(0.05,0.4,sprintf('%.2f < r < %.2f',r(ii*2-1),r(ii*2)),'FontSize',14);
text(0.05,0.2,'4 < q_{\psi} < 5','FontSize',14);
%text(0.035,0.2,sprintf('%.2f < \\rho < %.2f',r(ll),r(ll+1)),'FontSize',14);
if ii==3; ax.XTick = [0.05 0.1 0.2 0.5 1 2 5 10]; xlabel('\nu_{eff}');
else ax.XTickLabel = [];end

end
end
end

%% beta_BB vs. nu_eff in different r, Nl, Power with P_ICRH and P_LH
% for PoP
%r=[-0.75 -0.55 -0.35 -0.15 0.05 0.20 0.40 0.60];%r=[-0.15 0.05];
r=[-0.55 -0.45 -0.15 -0.05 0.45 0.55];
q=[4 5];
Nl=[2 8]*1e19;
P=[1.5 3.5];
%P=[0.5 1.5 2.5 3.5 10];

for jj=1:length(P)-1
for kk=1:length(Nl)-1

figure('Pos',[60    85   558   590],'Color','w');  
for ii=1:length(r)/2
 subplot(length(r)/2,1,ii)
%subplot(121) % LH
ind=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
  d.PLH>P(jj) & d.PLH<P(jj+1) & d.PECRH<0.1 & d.PICRH<0.1 &...
   d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
plot(d.nu_eff_l(ind),d.betaBB(ind),'r^');hold on

ind2=indOK&d.rhoc>r(ii*2-1)&d.rhoc<r(ii*2)&d.qpsi>q(1)&d.qpsi<q(2)&...
 d.PICRH>P(jj) & d.PICRH<P(jj+1) & d.PECRH<0.1 & d.PLH<0.1 &...
 d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
 plot(d.nu_eff_l(ind2),d.betaBB(ind2),'bo');

ax=gca;
 set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[0.05 10],'YLim',[0.5 2.5],...
   'TickDir','out','XScale','log')
if ii==3; xlabel('\nu_{eff}'); else ax.XTickLabel = []; end
ax.Position(1)=0.15;
ax.Position(4)=0.25;
ylabel('\beta_{BB}') 
lg=legend('LH','ICRH','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(0.05,0.8,txt{ii},'FontSize',18)
text(0.05,0.4,sprintf('%.2f < r < %.2f',r(ii*2-1),r(ii*2)),'FontSize',14);
text(0.05,0.2,'4 < q_{\psi} < 5','FontSize',14);
%text(0.035,0.2,sprintf('%.2f < \\rho < %.2f',r(ll),r(ll+1)),'FontSize',14);
if ii==3; ax.XTick = [0.05 0.1 0.2 0.5 1 2 5 10]; xlabel('\nu_{eff}');
else ax.XTickLabel = [];end

end
end
end