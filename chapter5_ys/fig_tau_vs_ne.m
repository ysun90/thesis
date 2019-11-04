close all
clear

%% load database
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
% Ohmic 
indOH = (d.PICRH + d.PECRH + d.PLH)<0.1;
indLOC = d.Nla<2.6*d.Ip*0.9 & indOH; 
indSOC = d.Nla>2.6*d.Ip*1.1 & indOH;
indLS = {indLOC,indSOC};

%% Plot hist of Ip
figure('Pos',[331    89   956   561])
histogram(d.Ip(indOH & d.Nla>0),'BinWidth',0.01)
grid on
ax = gca;
ax.FontSize = 20;
ax.XTick = 0.4:0.1:1.3;
xlabel('I_{P} [MA]')
ylabel('Counts')
title('Histogram of I_{P} in Ohmic discharges')

%% tau vs ne
Ip = 0.5:0.1:1.2;
dIp = 0.01;
XLim = {[0 4],[0 4],[.5 4.5],[.5 4.5],[1 5],[1 5],[1 5],[1.5 5.5]};
YLim = {[.05 .25],[.1 .3],[.1 .3],[.1 .3],[.1 .3],[.1 .3],[.1 .3],[.15 .35]};
nLS = [1.4 1.6 1.85 2.2 2.5 2.6 2.9 3.3];
figure('color','w','Position',[29   100  1100  600])
for i = 1:8
subplot(2,4,i)
box on
hold on
ind = d.Ip>(Ip(i)-dIp) & d.Ip<(Ip(i)+dIp) & d.Nla>0 & indOH & d.EquiTag==1 & d.rhoc>-1 & d.rhoc<1;
x = d.Nla(ind);
y = d.Tau(ind);
plot(x,y,'.');
dr = 0.2;
xc = dr-0:2*dr:6-dr;
xl = xc-dr;
xu = xc+dr;
Y = zeros(length(xc),1);
%Y2 = zeros(length(xc),1);
for j = 1:length(xc)
    indx = x>xl(j) & x<xu(j);                
    if ~isempty(y(indx))
        Y(j) = median(y(indx));
        %Y2(j) = mean(y(indx));
    else
        Y(j) = NaN;
        %Y2(j) = NaN;
    end
end
plot(xc,Y,'r-','LineWidth',2)
%plot(xc,Y2,'y-','LineWidth',2)
plot([nLS(i) nLS(i)],YLim{i},'k--','LineWidth',2)
ax = gca;
ax.LineWidth=1.5;
ax.FontSize=12;
ax.Position(3) = 0.19;
ax.Position(4) = 0.35;
ax.XLim = XLim{i};
ax.YLim = YLim{i};
t = text(XLim{i}(2),YLim{i}(1),sprintf('I_{P} = %2.1f [MA]',Ip(i)));
t.VerticalAlignment = 'bottom';
t.HorizontalAlignment = 'right';
t.FontSize = 12;
hold off
end
xlabel('Central line averaged density <n_{l,0}> [10^{19} m^{-3}]',...
  'FontSize',18);
ylabel('Energy confinement time \tau_{E} [s]','FontSize',18);
lgd = legend('Data points','Median value','Position of the Knee point');
lgd.Orientation = 'horizontal';
lgd.Position = [0.9132    0.8907    0.0782    0.0915];
lgd.FontSize = 12;