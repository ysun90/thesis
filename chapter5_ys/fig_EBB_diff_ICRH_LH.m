%% Turbulence database by Taylor model
clear
% Reduced database
load('DB.mat'); d=DB;
% Averaged linear density
d.Nla = d.Nl./(2*d.a); % [10^19 m-3]
% Normalized radius
load('rhoc.mat'); d.rhoc=r;
% Taylor model
fitT = load('T358874.mat');
fitG = load('GG358874.mat');
d.EBB = fitT.EBB'; d.WBB=fitT.SDBB';
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
clear DB reflec fitT r
%d.EBB=d.EBB.*d.ES0;

% Filter data
indOK = d.FVAL<0.002 & lg(d.AmpBB./d.Noise)>25 & d.EquiTag==1 & ...
    abs(d.muBB)<50 & d.Nla>0 & d.f_plateau<1e3;
% Ohmic 
indOH = (d.PICRH + d.PECRH + d.PLH)<0.1;
indLOC = d.Nla<2.6*d.Ip*0.9 & indOH; 
indSOC = d.Nla>2.6*d.Ip*1.1 & indOH;
indLS = {indLOC,indSOC};
% PICH
indICRH1 = d.PICRH>0.5 & d.PICRH<1.5 & d.PECRH<0.1 & d.PLH<0.5;
indICRH2 = d.PICRH>1.5 & d.PICRH<2.5 & d.PECRH<0.1 & d.PLH<0.5;
indICRH3 = d.PICRH>2.5 & d.PICRH<3.5 & d.PECRH<0.1 & d.PLH<0.5;
indICRH4 = d.PICRH>3.5 & d.PICRH<inf & d.PECRH<0.1 & d.PLH<0.5;
indICRH = {indICRH1,indICRH2,indICRH3,indICRH4};
% PLH
indLH1 = d.PLH>0.5 & d.PLH<1.5 & d.PECRH<0.1 & d.PICRH<0.5;
indLH2 = d.PLH>1.5 & d.PLH<2.5 & d.PECRH<0.1 & d.PICRH<0.5;
indLH3 = d.PLH>2.5 & d.PLH<3.5 & d.PECRH<0.1 & d.PICRH<0.5;
indLH4 = d.PLH>3.5 & d.PLH<inf & d.PECRH<0.1 & d.PICRH<0.5;
indLH = {indLH1,indLH2,indLH3,indLH4};
% PECRH
indECRH=d.PLH<0.5 & d.PECRH>0.1 & d.PICRH<0.5;
indECRH = {indECRH,indECRH,indECRH,indECRH};
% Safety factor 
q = [3 4 5 6 10];
indq = cell(1,length(q)-1);
for i = 1:length(indq)
    indq{i} = d.qpsi>q(i) & d.qpsi<q(i+1);
end
% Subinterval of radil positions
dr = 0.05;rc = dr-1:2*dr:1-dr;
rl = rc-dr;rr = rc+dr;
inddr = cell(1,length(rc));
for i = 1:length(inddr)
    inddr{i} = d.rhoc>=rl(i) & d.rhoc<rr(i);
end
%VAR = zeros(1,length(rc));
R2 = zeros(1,length(rc));

%% show the remarkable difference of E_BB with ICRH and LH
% with histogram

%%% PARAMETERS %%%
qmin = 3; qmax = 10; % edge safety factor
rin = 0; rout = 0.4; dr = 0.05; % radial position in- and outside q=1

indq = d.qpsi > qmin & d.qpsi < qmax;

%%% PLOT %%%
f1 = figure; set(f1,'color','w','position',[158 196 982 416])

%%% RATIO OF E_BB INSIDE AND OUTSIDE q=1 %%%
s1 = subplot(121);

indrin = abs(d.rhoc-rin)<=dr;
indrout = abs(d.rhoc-rout)<=dr;

for ii=1:length(indICRH)
  
  EBBinICRH(ii) = median(d.EBB(indq & indOK & indrin & indICRH{ii}));
  EBBoutICRH(ii) = median(d.EBB(indq & indOK & indrout & indICRH{ii}));
  EBBinLH(ii) = median(d.EBB(indq & indOK & indrin & indLH{ii}));
  EBBoutLH(ii) = median(d.EBB(indq & indOK & indrout & indLH{ii}));
  
end

hold on
p1 = plot(EBBinICRH./EBBoutICRH,'bs--','MarkerSize',18,'MarkerFaceColor','b');
p2 = plot(EBBinLH./EBBoutLH,'r^--','MarkerSize',18,'MarkerFaceColor','r');
hold off

xlabel('Heat. Power [MW]'); ylabel('E_{BB}^{core}/E_{BB}^{edge}'); 
set(s1,'LineWidth',1.5,'FontSize',18,'XLim',[0 5],'YLim',[0 1.2],'box','on')
legend([p1,p2],'ICRH','LH');
text(1,1,'(a)','FontSize',16);

%%% HISTOGRAM OF CONFINEMENT WITH ICRH AND LH
s2 = subplot(122);

indicrh = d.PICRH>.1 & d.PECRH<.1 & d.PLH<.1;
indlh = d.PLH>.1 & d.PECRH<.1 & d.PICRH<.1;
indq = d.qpsi > qmin & d.qpsi < qmax;

hold on
h1 = histogram(d.Tau(indq & indOK & indicrh),'BinWidth',.01,'Normalization','probability');
h2 = histogram(d.Tau(indq & indOK & indlh),'BinWidth',.01,'Normalization','probability');
hold off

xlabel('\tau_{E} [s]'); ylabel('Number of spectra'); 
set(s2,'LineWidth',1.5,'FontSize',18,'XLim',[0 .3],'box','on')
legend([h1,h2],'ICRH','LH');
text(0,3000,'(b)','FontSize',16);
%text(0,3000,sprintf('%d < q_{\\psi} < %d',qmin,qmax),'FontSize',14);


%% show the remarkable difference of E_BB with ICRH and LH

%%% PARAMETERS %%%
qmin = 3; qmax = 8; % edge safety factor
rin = 0; rout = 0.4; dr = 0.05; % radial position in- and outside q=1

indq = d.qpsi > qmin & d.qpsi < qmax;

%%% PLOT %%%
f1 = figure; set(f1,'color','w','position',[158 196 982 416])

%%% RATIO OF E_BB INSIDE AND OUTSIDE q=1 %%%
s1 = subplot(121);

indrin = abs(d.rhoc-rin)<=dr;
indrout = abs(d.rhoc-rout)<=dr;

for ii=1:length(indICRH)
  
  EBBinICRH(ii) = median(d.EBB(indq & indOK & indrin & indICRH{ii}));
  EBBoutICRH(ii) = median(d.EBB(indq & indOK & indrout & indICRH{ii}));
  EBBinLH(ii) = median(d.EBB(indq & indOK & indrin & indLH{ii}));
  EBBoutLH(ii) = median(d.EBB(indq & indOK & indrout & indLH{ii}));
  
end

% hold on
% p1 = plot(EBBinICRH./EBBoutICRH,'bs--','MarkerSize',18,'MarkerFaceColor','b');
% p2 = plot(EBBinLH./EBBoutLH,'r^--','MarkerSize',18,'MarkerFaceColor','r');
% hold off

hold on
p1 = plot(EBBinICRH,'bs--','MarkerSize',18,'MarkerFaceColor','b');
p2 = plot(EBBinLH,'r^--','MarkerSize',18,'MarkerFaceColor','r');
hold off

xlabel('Heat. Power [MW]'); ylabel('E_{BB}^{core}'); 
set(s1,'LineWidth',1.5,'FontSize',18,'XLim',[0 5],'YLim',[0 1.2],'box','on')
legend([p1,p2],'ICRH','LH');
text(1,1,'(a)','FontSize',16);

%%% HISTOGRAM OF CONFINEMENT WITH ICRH AND LH
s2 = subplot(122);

indicrh = d.PICRH>.1 & d.PECRH<.1 & d.PLH<.1;
indlh = d.PLH>.1 & d.PECRH<.1 & d.PICRH<.1;
indq = d.qpsi > qmin & d.qpsi < qmax;

for ii=1:length(indICRH)
  
  tauICRH(ii) = median(d.Tau(indq & indOK & indICRH{ii}));
  tauLH(ii) = median(d.Tau(indq & indOK & indLH{ii}));
  
end

hold on
p1 = plot(tauICRH,'bs--','MarkerSize',18,'MarkerFaceColor','b');
p2 = plot(tauLH,'r^--','MarkerSize',18,'MarkerFaceColor','r');
hold off

xlabel('Heat. Power [MW]'); ylabel('\tau_E [s]'); 
set(s2,'LineWidth',1.5,'FontSize',18,'XLim',[0 5],'YLim',[0 .2],'box','on')
legend([p1,p2],'ICRH','LH');
text(1,.1,'(b)','FontSize',16);
%text(0,3000,sprintf('%d < q_{\\psi} < %d',qmin,qmax),'FontSize',14);


