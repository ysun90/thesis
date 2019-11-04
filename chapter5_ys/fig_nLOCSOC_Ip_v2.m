%% n vs Ip in TS database
load('DB.mat'); d=DB;
d.Nla = d.Nl./(2*d.a); % [10^19 m-3]
indOH = (d.PICRH + d.PECRH + d.PLH)<0.1;
% Fit nLOCSOC vs Ip in TS
Ip = [0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2]; % MA
IpNeg = ones(1,7)*0.1;
IpPos = IpNeg;
nLOCSOC = [1.4 1.6 1.85 2.2 2.5 2.6 2.9 3.2]; % 10^19 [m^-3]
%nLOCSOC = [1.5 1.7 1.85 2.1 2.5 2.6 3.2]; % 10^19 [m^-3]
nNeg = [0.1 0.2 0.1 0.2 0.2 0.1 0.1 0.1]; % uncertainties
nPos = nNeg;

figure('Position',[169 71 1105 562],'Color','w');
Ip = 0.5:0.1:1.2;
dIp = 0.01;
XLim = {[0 4],[0 4],[.5 4.5],[.5 4.5],[1 5],[1 5],[1 5],[1.5 5.5]};
YLim = {[.05 .25],[.1 .3],[.1 .3],[.1 .3],[.1 .3],[.1 .3],[.1 .3],[.15 .35]};
nLS = [1.4 1.6 1.85 2.2 2.5 2.6 2.9 3.2];
for i = 1:8
  subplot(2,4,i)
box on
hold on
ind = d.Ip>(Ip(i)-dIp) & d.Ip<(Ip(i)+dIp) & d.Nla>0 & indOH & d.EquiTag==1 & d.rhoc>-1 & d.rhoc<1;
kappa(i)=mean(d.Ellip(ind));
a(i)=mean(d.a(ind));
x = d.Nla(ind);
y = d.Tau(ind);
plot(x,y,'.');
dr = 0.2;
xc = dr-0:2*dr:6-dr;
xl = xc-dr;
xu = xc+dr;
Y = zeros(length(xc),1);
Y2 = zeros(length(xc),1);
for j = 1:length(xc)
    indx = x>xl(j) & x<xu(j);                
    if ~isempty(y(indx))
        Y(j) = median(y(indx));
        Y2(j) = mean(y(indx));
    else
        Y(j) = NaN;
        Y2(j) = NaN;
    end
end
plot(xc(1:end-2),Y(1:end-2),'r-','LineWidth',2)
%plot(xc,Y2,'y-','LineWidth',2)
plot([nLS(i) nLS(i)],YLim{i},'k--','LineWidth',2)
ax = gca;
ax.FontSize=14;
ax.LineWidth=1.5;
ax.XLim = XLim{i};
ax.YLim = YLim{i};
t = text(XLim{i}(2),YLim{i}(1),sprintf('I_{P} = %2.1f [MA]',Ip(i)));
t.FontSize=14;
hold off
end
xlabel('n_{e} [10^{19} m^{-3}]');
ylabel('\tau_{E} [s]');
lgd = legend('Data points','Median value','n_{LOC-SOC}');
lgd.FontSize = 12;
text(1,.25,'(a)','FontSize',16)

%% Fit by nLOCSOC = 2.6*Ip
I = 0.5:0.1:1.2;
y = 2.6*I;
figure('Position',[169 71 1105 562],'Color','w')
hold on;
grid on;
errorbar(Ip,nLOCSOC,nPos,'ko','MarkerSize',10,'MarkerFaceColor','k')
l2=plot(Ip,y,'r-','LineWidth',3);
p1=plot(0.5,1.5,'mv','MarkerSize',10,'MarkerFaceColor','r');
p2=plot(1.2,3,'m^','MarkerSize',10,'MarkerFaceColor','r');
ax = gca;
ax.Box = 'on';
ax.LineWidth=1.5;
ax.FontSize = 14;
ax.XLim = [0.4 1.3];
ax.YLim = [1 3.5];
xlabel('I_{p} [MA]');
ylabel('n_{LOC-SOC} [10^{19} m^{-3}]');
%title('LOC-SOC Density Threshold');
lg = legend([l2,p1,p2],{'Linear fit ($R^2=0.98$)',...
  '#41261 - 41272 (Ip = 0.5 MA)','#41003 - 41013 (Ip = 1.2 MA)'},...
  'Interpreter','latex');
lg.FontSize = 12;
lg.Location = 'best';
text(.5,3,'(b)','FontSize',16)


% %%
% % Fit nLOCSOC vs Ip in TS
% IpNeg = ones(1,7)*0.1;
% IpPos = IpNeg;
% nLOCSOC = [1.4 1.6 1.85 2.1 2.5 2.6 2.9 3.2]; % 10^19 [m^-3]
% %nLOCSOC = [1.5 1.7 1.85 2.1 2.5 2.6 3.2]; % 10^19 [m^-3]
% nNeg = [0.1 0.2 0.1 0.2 0.2 0.1 0.1 0.1]; % uncertainties
% nPos = nNeg;

% % Fit by nLOCSOC = 2.6*Ip
% I = 0.5:0.1:1.2;
% x=I.*sqrt(kappa)./a.^2;
% y = 1.33*x;
% figure('Position',[169 71 1105 562],'Color','w');
% hold on;
% grid on;
% errorbar(x,nLOCSOC,nPos,'k.','MarkerSize',20)
% l2=plot(x,y,'r-','LineWidth',3);
% ax = gca;
% ax.Box = 'on';
% ax.FontSize = 20;
% ax.XLim = [0.9 2.4];
% ax.YLim = [1 3.5];
% xlabel('$I_{p}\sqrt\kappa/a^2$ [MA $\cdot$ $m^{-2}]$','Interpreter','latex');
% ylabel('n_{LOC-SOC} [10^{19} m^{-3}]');
% %title('LOC-SOC Density Threshold');
% lg = legend(l2,{'Fit ($R^2=0.98$)'},'Interpreter','latex');
% lg.FontSize = 20;
% lg.Location = 'best';