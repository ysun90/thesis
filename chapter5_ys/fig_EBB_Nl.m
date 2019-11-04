%% Turbulence database by Taylor model

% Reduced database
load('DB.mat'); d=DB;
db = load('Database_v10_plasma_temp.mat');
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
d.nu_eff_l=0.1*2*d.Rp.*(d.Nec/1e19)./(d.Tec.^2);
%d.nu_eff_l=d.Rp./(d.a./d.Rp).*(d.Nec./1e19)./(d.Tec.^2);
%d.nu_eff_l=d.Rp./(0.1./d.Rp).*(d.Nec./1e19)./(d.Tec.^2);
%d.nu_eff_l=0.0468*d.Rp.*(d.Nec/1e19)./(d.Tec.^2);
clear DB reflec fitT r Nelc Nelc_robust db
%d.EBB=d.EBB.*d.ES0;

% Filter data
 %indOK = d.FVAL<0.002 & lg(d.EBB./(d.Noise*1025))>0 & d.EquiTag==1 & ...
    % abs(d.muBB)<50;
  %indOK = d.FVAL<0.002 & lg(d.EBB./(d.Noise*1025))>0 & d.EquiTag==1 & ...
   % abs(d.muBB)<50 & d.Tec>0;
 indOK = d.FVAL<0.002 & lg(d.AmpBB./d.Noise)>25 & d.EquiTag==1 & ...
    abs(d.muBB)<50 & d.Nla>0 & d.f_plateau<1e3;
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
txt={'(a)','(b)','(c)'};

%% E_BB vs. nl in different q, Nl, Power with P_ICRH and P_LH
% for PoP
%r=[-0.75 -0.55 -0.35 -0.15 0.05 0.20 0.40 0.60];%r=[-0.15 0.05];
r=[-0.15 0.05];
q=[3 4 5 6];
Nl=[1 10]*1e19;
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
plot(d.Nl(ind),d.EBB(ind),'r^');hold on

%ax=gca;
% set(ax,'FontSize',14,'LineWidth',1.5,'XLim',[0.03 2],'YLim',[0 1],'XScale','log')
%title(sprintf('%d < q < %d',q(ii),q(ii+1)))
% title(sprintf('%d < q < %d, %.2g < Nl < %.2g, %.2g < P < %.2g, %.2g < r/a < %.2g',...
%   q(ii),q(ii+1),Nl(kk),Nl(kk+1),P(jj),P(jj+1),r(1),r(2)))

ind2=indOK&d.rhoc>r(ll)&d.rhoc<r(ll+1)&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
 d.PICRH>P(jj) & d.PICRH<P(jj+1) & d.PECRH<0.1 & d.PLH<0.1 &...
 d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
 plot(d.Nl(ind2),d.EBB(ind2),'bo');

ax=gca;
 set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[2 8],'YLim',[0 1],...
   'TickDir','out','XScale','linear')
ax.Position(1)=0.15;
ax.Position(4)=0.24;
ylabel('E_{BB}') 
lg=legend('LH','ICRH','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(3,0.8,txt{ii},'FontSize',18)
text(5,0.4,sprintf('%d < q_{\\psi} < %d',q(ii),q(ii+1)),'FontSize',14);
text(8,0.2,'\rho = - 0.1 \pm 0.05','FontSize',14);
%text(0.035,0.2,sprintf('%.2f < \\rho < %.2f',r(ll),r(ll+1)),'FontSize',14);
  if ii==3
    xlabel('N_l [\times 10^{19} m^{-3}]'); %ax.XTickLabel = {'0.1','1'};
  else
    ax.XTickLabel = [];
  end
end
end
end
end

%% E_BB vs. Nl in different q, Nl, with OH
% for PoP
r=[-0.15 -0.05];
q=[3 4 5 6];
Nl=[0 10]*1e19;
%P=[1.5 3.5];
%P=[0.5 1.5 2.5 3.5 10];

for kk=1:length(Nl)-1
figure('Pos',[60    85   558   590],'Color','w');  
for ii=1:length(q)-1
subplot(length(q)-1,1,ii)
ind=indOK&d.rhoc>r(1)&d.rhoc<r(2)&d.qpsi>q(ii)&d.qpsi<q(ii+1)&...
  d.PLH<0.1 & d.PECRH<0.1 & d.PICRH<0.1 &...
   d.Nelc>Nl(kk) & d.Nelc<Nl(kk+1);
hold on
plot(d.Nl(ind&indLOC),d.EBB(ind&indLOC),'r^');
plot(d.Nl(ind & ~indLOC & ~indSOC),d.EBB(ind & ~indLOC & ~indSOC),'kd');
plot(d.Nl(ind&indSOC),d.EBB(ind&indSOC),'bo');
hold off

ax=gca;
set(ax,'FontSize',18,'LineWidth',1.5,'XLim',[1 8],'YLim',[0 1],...
   'TickDir','out','XScale','linear')
box on
%grid on
ax.Position(1)=0.15;
ax.Position(4)=0.25;
ylabel('E_{BB}') 
lg=legend('LOC','Transition','SOC','Location','southeast');
lg.FontSize=14;
lg.Orientation = 'horizontal';
text(0.035,0.8,txt{ii},'FontSize',18)
text(0.035,0.4,sprintf('%d < q_{\\psi} < %d',q(ii),q(ii+1)),'FontSize',14);
text(0.035,0.2,'\rho = - 0.1 \pm 0.05','FontSize',14);
  if ii==3
    xlabel('N_l [\times 10^{19} m^{-3}]'); %ax.XTickLabel = {'0.1','1'};
  else
    ax.XTickLabel = [];
  end
end
end
% for i=1:48; figure(i); saveas(gcf,['f',num2str(i)],'png'); end