
%% Turbulence database by Taylor model

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


%% ICRH
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
    figure('Position',[58          60        1253         616],'Color','w');
    for i = 1:m
        for j = 1:n
            ind = indOK & indq{i} & indICRH{j};
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
            subplot(m,n,(n*(i-1)+j));
            hold on
            % Data points
            plot(d.rhoc(ind),d.(VarNames{k})(ind),'.','Color',[0 1 1]);
            % Median values
            ax2 = plot(rc,R.(VarNames{k})(:,sub2ind([m,n],i,j)),'rs--','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');
            % q=1 surface
            rq1 = mean(1./d.qpsi(ind));
            rq1_md(i) = rq1;
            l1 = line([-rq1 -rq1],[YLimICRH{k}]);
            l1.Color = 'k';
            l1.LineWidth = 2;
            l1.LineStyle = '--';
            l2 = line([rq1 rq1],[YLimICRH{k}]);
            l2.Color = 'k';
            l2.LineWidth = 2;
            l2.LineStyle = '--';
            % LOC
            %ax3=plot(rc,VAR.(VarNames{k})(:,i),'r-','LineWidth',2);
            % SOC
            %ax4=plot(rc,VAR.(VarNames{k})(:,i+4),'k:','LineWidth',2);
            % OH
            ax5=plot(rc,(VAR.(VarNames{k})(:,i)+VAR.(VarNames{k})(:,i+4))/2,'ko--','LineWidth',2);
            hold off
            
            ax = gca;
            ax.Position(3) = 0.19;
            ax.Position(4) = 0.19;
            ax.Box = 'on';
            ax.XLim = [-1 0.6];
            ax.XTick = -1:0.5:1;
            ax.YLim = YLimICRH{k};
            ax.TickDir = 'out';
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
            % YTICK
            if j~=1 
                ax.YTickLabel = [];
            end        
            % TITLE and YLABEL
            if i==2
                title(tP{j},'FontSize',14);
                if j==1
                    ylabel(yl{k},'FontSize',14);
                end
            end
            
        end
    end
end

%% LH
VarNames = {'EBB','k2DBB','tauBB','kutauBB','muBB','ELF','sigmaLF','muLF','ECS','sigmaCS','ES0'};
tq = {'2 < q_{\psi} < 4','4 < q_{\psi} < 5','5 < q_{\psi} < 6','6 < q_{\psi} < 8'};
tP = {'0.5 < P_{LH} < 1.5 (MW)','1.5 < P_{LH} < 2.5 (MW)','2.5 < P_{LH} < 3.5 (MW)','P_{LH} > 3.5 (MW)'};
YLimICRH = {[0 1],[0 10],[0,80],[0 10],[-50 50],[0 1],[0 20],[-5 5],[0 .6],[0 1.5],[0 60]};
yl = {'E_{BB}','k2D_{BB}','\tau_{BB}','ku\tau','\mu_{BB}',...
    'E_{LF}','\sigma_{LF}','\mu_{LF}','E_{CS}','\sigma_{CS}','E_{S0}',};
l = length(VarNames);
m = length(indq);
n = length(indLH);
rq1_md = zeros(1,m);
hold on
for k = 1:1
    %figure('Position',[58          60        1253         616],'Color','w');
    for i = 2
        for j = 1:n
            ind = indOK & indq{i} & indLH{j};
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
            subplot(m,n,(n*(i-1)+j));
            hold on
            % Data points
            plot(d.rhoc(ind),d.(VarNames{k})(ind),'.','Color',[0 1 1]);
            % Median values
            ax6 = plot(rc,R.(VarNames{k})(:,sub2ind([m,n],i,j)),'b*--','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','g');
            % q=1 surface
            rq1 = mean(1./d.qpsi(ind));
            rq1_md(i) = rq1;
            l1 = line([-rq1 -rq1],[YLimICRH{k}]);
            l1.Color = 'k';
            l1.LineWidth = 2;
            l1.LineStyle = '--';
            l2 = line([rq1 rq1],[YLimICRH{k}]);
            l2.Color = 'k';
            l2.LineWidth = 2;
            l2.LineStyle = '--';
%             % LOC
%             ax3=plot(rc,VAR.(VarNames{k})(:,i),'r-','LineWidth',2);
%             % SOC
%             ax4=plot(rc,VAR.(VarNames{k})(:,i+4),'k:','LineWidth',2);
            ax5=plot(rc,(VAR.(VarNames{k})(:,i)+VAR.(VarNames{k})(:,i+4))/2,'ko--','LineWidth',2);
            hold off
            
            ax = gca;
            ax.Position(3) = 0.19;
            ax.Position(4) = 0.19;
            ax.Box = 'on';
            ax.XLim = [-1 0.6];
            ax.XTick = -1:0.5:1;
            ax.YLim = YLimICRH{k};
            ax.TickDir = 'out';
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
            % YTICK
            if j~=1 
                ax.YTickLabel = [];
            end        
            % TITLE and YLABEL
            if i==2
                title(tP{j},'FontSize',14);
                if j==1
                    ylabel(yl{k},'FontSize',14);
                end
            end
            
        end
    end
            legend([ax5 ax2 ax6 l1],'Ohmic','ICRH','LH','q=1 surface');
end