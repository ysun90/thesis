
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
d.EBB = fitT.EBB'; d.EBBc = d.EBB./sqrt(d.Lepsi.*d.f_plateau);
d.EBB = d.EBBc;
d.WBB=fitT.SDBB';
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

% Filter data
indOK = d.FVAL<0.002 & lg(d.AmpBB./d.Noise)>25 & d.EquiTag==1 & ...
    abs(d.muBB)<50 & d.Nla>0 & d.f_plateau<1e3;
% Ohmic 
indOH = (d.PICRH + d.PECRH + d.PLH)<0.1;
indLOC = d.Nla<2.6*d.Ip*0.9 & indOH; 
indSOC = d.Nla>2.6*d.Ip*1.1 & indOH;
indLS = {indLOC,indSOC};
% PICH
indICRH = d.PICRH>4 & d.PECRH<0.1 & d.PLH<0.5;
% PLH
indLH = d.PLH>3 & d.PECRH<0.1 & d.PICRH<0.5;
% Mix ICRH and LH
indMix= d.PICRH>0.5 & d.PICRH<1.5 & d.PLH>2.5 & d.PLH<3.5 & d.PECRH<0.1 ;
% index all
indALL={indOH,indICRH,indLH};
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

%% Mix ICRH and LH
VarNames = {'EBB','WBB','k2DBB','tauBB','kutauBB','muBB','ELF','sigmaLF','muLF','ECS','sigmaCS','ES0','Lepsi','betaBB'};
tq = {'2 < q_{\psi} < 4','4 < q_{\psi} < 5','5 < q_{\psi} < 6','6 < q_{\psi} < 8'};
tP = {'0.5 < P_{ICRH} < 1.5 (MW)','1.5 < P_{ICRH} < 2.5 (MW)','2.5 < P_{ICRH} < 3.5 (MW)','P_{ICRH} > 3.5 (MW)'};
YLimICRH = {[0 1],[0 120],[0 10],[0,80],[0 10],[-50 50],[0 1],[0 20],[-5 5],[0 .6],[0 1.5],[0 60],[0 0.5],[.5 2.5]};
yl = {'E_{BB}','k2D_{BB}','\tau_{BB}','ku\tau','\mu_{BB}',...
    'E_{LF}','\sigma_{LF}','\mu_{LF}','E_{CS}','\sigma_{CS}','E_{S0}'};
l = length(VarNames);
m = length(indq);
n = length(indICRH);
rq1_md = zeros(1,m);
for k = 1:1
    %figure('Position',[58          60        1253         616],'Color','w');
    for i = 2
        for j = 1:3
            ind = indOK & indq{i} & indALL{j};
            numSHOT=numel(unique(d.Shot(ind)));
            numSpec=length(find(ind));
            % Averaged values vs position
            for ii = 1:length(rc)
                indr = d.rhoc>rl(ii) & d.rhoc<rr(ii);
                indR = ind & indr;
                R_temp = d.(VarNames{k})(indR);
                if ~isempty(R_temp)
                    R.(VarNames{k})(ii,sub2ind([m,n],i,j)) = median(R_temp);
                else
                    R.(VarNames{k})(ii,sub2ind([m,n],i,j)) = NaN;
                end
            end
            %subplot(m,n,(n*(i-1)+j));
            hold on
            % Data points
            %plot(d.rhoc(ind),d.(VarNames{k})(ind),'.','Color',[0 1 1]);
            % Median values
            %ax2 = plot(rc,R.(VarNames{k})(:,sub2ind([m,n],i,j)),'rs--','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');
            % q=1 surface
            rq1 = mean(1./d.qpsi(ind));
            rq1_md(j) = rq1;
            nshot(j)=numSHOT;
            nspec(j)=numSpec;
%             l1 = line([-rq1 -rq1],[YLimICRH{k}]);
%             l1.Color = 'k';
%             l1.LineWidth = 2;
%             l1.LineStyle = '--';
%             l2 = line([rq1 rq1],[YLimICRH{k}]);
%             l2.Color = 'k';
%             l2.LineWidth = 2;
%             l2.LineStyle = '--';
            % LOC
            %ax3=plot(rc,VAR.(VarNames{k})(:,i),'r-','LineWidth',2);
            % SOC
            %ax4=plot(rc,VAR.(VarNames{k})(:,i+4),'k:','LineWidth',2);
            % OH
            %ax5=plot(rc,(VAR.(VarNames{k})(:,i)+VAR.(VarNames{k})(:,i+4))/2,'ko--','LineWidth',2);
            hold off
            
            ax = gca;
%             ax.Position(3) = 0.19;
%             ax.Position(4) = 0.19;
            ax.Box = 'on';
            ax.XLim = [-1 0.6];
            ax.XTick = -1:0.5:1;
            ax.YLim = YLimICRH{k};
            ax.FontSize=14;
            ax.LineWidth=1.5;
            ql = [3,4,5,6];
            qu = [4,5,6,10];
            if j==n
                t = text(1,YLimICRH{k}(1)*2,sprintf('%d < q_{\\psi} < %d',ql(i),qu(i)));
                t.FontSize = 14;
                t.FontWeight = 'Bold';
                t.Margin = 0.01;
                t.Rotation = -90;
                %t.BackgroundColor = 'white';
                t.VerticalAlignment = 'top';
                t.HorizontalAlignment = 'right'; 
            end
            % XTICK
            if i==2
                ax.XTickLabel = {'-1','-0.5','0','0.5'};
                xlabel('\rho');
            else
                ax.XTickLabel = [];
            end

            
        end
    end
            ylabel('E_{BB}')
end

%% Plot;
h=figure;
set(h,'Color','white')
EBB=R.EBB(:,[2,6,10]);
ls={'ko--','rs--','g^--','b*--'};
for ii=1:3
    hold on
    plot(rc(2:end),EBB(2:end,ii),ls{ii},'LineWidth',2,'MarkerSize',16)
    hold off   
end
set(gca,'Box','on','LineWidth',2,'FontSize',20,'XLim',[-0.9 0.5],'YLim',[0 0.2])
xlabel('\rho');
ylabel('corrected E_{BB}')
legend({'Ohmic - 29463 spectra in 1384 shots',...
       'P_{ICRH} > 4 MW  - 851 spectra in 55 shots',...
       'P_{LH} > 3 MW - 855 spectra in 35 shots'},'FontSize',16)
hold on
plot([-rq1_md(1) -rq1_md(1)],[0 0.05],'k-.','LineWidth',2)
plot([rq1_md(1) rq1_md(1)],[0 0.05],'k-.','LineWidth',2)
plot([-rq1_md(2) -rq1_md(2)],[0.1 0.2],'r-.','LineWidth',2)
plot([rq1_md(2) rq1_md(2)],[0.1 0.2],'r-.','LineWidth',2)
plot([-rq1_md(3) -rq1_md(3)],[0.05 0.1],'g-.','LineWidth',2)
plot([rq1_md(3) rq1_md(3)],[0.05 0.1],'g-.','LineWidth',2)
% plot([-rq1_md(4) -rq1_md(4)],[0.7 0.9],'b-.','LineWidth',2)
% plot([rq1_md(4) rq1_md(4)],[0.7 0.9],'b-.','LineWidth',2)
hold off

text(-0.2,0.5,'q=1','Rotation',90,'FontSize',16,'Color','k')
text(-0.2,0.5,'q=1','Rotation',90,'FontSize',16,'Color','g')
text(-0.2,0.5,'q=1','Rotation',90,'FontSize',16,'Color','r')
text(0.2,0.5,'q=1','Rotation',90,'FontSize',16,'Color','k')
text(0.2,0.5,'q=1','Rotation',90,'FontSize',16,'Color','g')
text(0.2,0.5,'q=1','Rotation',90,'FontSize',16,'Color','r')
% text(0,0.5,'q=1','Rotation',90,'FontSize',16,'Color','b')