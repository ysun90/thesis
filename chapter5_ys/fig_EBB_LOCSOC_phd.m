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
fitG=load('GG358874.mat');
d.EBB = fitT.EBB'; d.WBB=fitT.SDBB';
d.betaBB=fitG.betaBB';
d.k2DBB=fitT.k2DBB'; %(k^2*u^2*tau)
d.tauBB=fitT.tauBB';d.kutauBB=sqrt(d.k2DBB.*d.tauBB);
d.ELF = fitT.ELF';d.sigmaLF=fitT.sigmaLF';d.muLF = fitT.muLF';
d.ECS=fitT.ECS';d.sigmaCS=fitT.sigmaCS';
d.ES0=fitT.ES0';
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
indICRH3 = d.PICRH>2.5 & d.PICRH<4 & d.PECRH<0.1 & d.PLH<0.5;
indICRH4 = d.PICRH>4 & d.PICRH<inf & d.PECRH<0.1 & d.PLH<0.5;
indICRH = {indICRH1,indICRH2,indICRH3,indICRH4};
% PLH
indLH1 = d.PLH>0.5 & d.PLH<1.5 & d.PECRH<0.1 & d.PICRH<0.1;
indLH2 = d.PLH>1.5 & d.PLH<2.5 & d.PECRH<0.1 & d.PICRH<0.1;
indLH3 = d.PLH>2.5 & d.PLH<3.5 & d.PECRH<0.1 & d.PICRH<0.1;
indLH4 = d.PLH>3.5 & d.PLH<inf & d.PECRH<0.1 & d.PICRH<0.1;
indLH5 = d.PLH>0.5 & d.PLH<inf & d.PECRH<0.1 & d.PICRH<0.1;
indLH = {indLH1,indLH2,indLH3,indLH4,indLH5};
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
%indF={'EBB','WBB','kutauBB','betaBB'};

%% LOC & SOC original profiles
% Names, titles, limits, labels in loop
VarNames = {'EBB','WBB','k2DBB','tauBB','kutauBB','muBB','ELF','sigmaLF','muLF','ECS','sigmaCS','ES0','Lepsi','betaBB','muCS'};
YLim = {[0 1],[0 120],[0 10],[0,80],[0 10],[-50 50],[0 1],[0 20],[-5 5],[0 .6],[0 1.5],[0 60],[0 0.5],[0.5 2.5],[-2 2]};
tVar = {'E_{BB}','k2D_{BB}','\tau_{BB}','ku\tau','\mu_{BB}',...
    'E_{LF}','\sigma_{LF}','\mu_{LF}','E_{CS}','\sigma_{CS}','E_{S0}','L_{\epsilon}','\beta_{BB}'};
tLS = {'LOC','SOC'};
l = length(VarNames);m = length(indq);n = length(indLS);
rq1_md = zeros(1,m);
for k = 1:1
        figure('Position',[480    49   841   635],'Color','w');
        for i = 1:m
            for j = 1:n
            ind = indOK & indq{i} & indLS{j};   
            
            % Median values within subintervals
            for ii = 1:length(rc)                
                VARdr = d.(VarNames{k})(ind & inddr{ii});
                if ~isempty(VARdr)
                    VAR.(VarNames{k})(ii,sub2ind([m,n],i,j)) = median(VARdr);
                else
                    VAR.(VarNames{k})(ii,sub2ind([m,n],i,j)) = NaN;
                end
            end 
            
            ax1 = subplot(m,n,(n*(i-1)+j));
            if j==1
                ax1.Position(1)=0.05;
            elseif j==2
                ax1.Position(1)=0.52;
            end
            ax1.Position(3) = 0.43;
            ax1.Position(4) = 0.19;
            hold on
            %plot(d.rhoc(ind),d.(VarNames{k})(ind),'.','Color',[0 1 1]);
            if j==1
                plot(d.rhoc(ind),d.(VarNames{k})(ind),'r.');
            elseif j==2
               plot(d.rhoc(ind),d.(VarNames{k})(ind),'b.');
            end
            
            rq1 = mean(1./d.qpsi(ind));
            rq1_md(i) = rq1;
            l1 = line([-rq1 -rq1],[YLim{k}]);
            l1.Color = 'k';
            l1.LineWidth = 2;
            l1.LineStyle = '--';
            l2 = line([rq1 rq1],[YLim{k}]);
            l2.Color = 'k';
            l2.LineWidth = 2;
            l2.LineStyle = '--';

            hold off
            
            ax = gca;
            ax.LineWidth=2;
            ax.Box = 'on';
            ax.XLim = [-1 0.6];
            ax.TickDir = 'out';
            ax.YLim = YLim{k};
            ax.FontSize = 14;
            
            if j==n; ax.YAxisLocation='right'; end
            
            if i==m
                ax.XTickLabel = {'-1','-0.5','0','0.5','1'};
                xlabel('\rho','FontSize',16);              
            else
                ax.XTickLabel = [];
            end
            
            if i==round(m/2) && j==1;  ylabel(tVar{k},'FontSize',16); end
            
           if i==1; title([tLS{j}],'FontSize',16); end
                     
            t = text(-0.95,YLim{k}(2)*0.95,sprintf('%d < q_{\\psi} < %d',q(i),q(i+1)));
            t.FontSize = 12;
            t.BackgroundColor = 'white';
            t.VerticalAlignment = 'top';
            t.HorizontalAlignment = 'left';  
            end
        end
%         lgd = legend(ax2,'Median values');
%         lgd.Position = [0.4382    0.0073    0.1950    0.0463];
%         lgd.Orientation = 'horizontal';
%         lgd.FontSize = 14;
end

 ylabel('E_{BB}');
% title('LOC')
% title('SOC')

%% LOC & SOC
% Names, titles, limits, labels in loop
txt={'(a)','(b)','(c)','(d)'};
VarNames = {'EBB','WBB','k2DBB','tauBB','kutauBB','muBB','ELF','sigmaLF','muLF','ECS','sigmaCS','ES0','Lepsi','betaBB'};
YLim = {[0 1],[0 120],[0 10],[0,80],[0 10],[-50 50],[0 1],[0 20],[-5 5],[0 .7],[0 1.5],[0 60],[0 0.5],[0.5 2.5]};
tVar = {'E_{BB}','w_{BB}','k2D_{BB}','\tau_{BB}','ku\tau','\mu_{BB}',...
    'E_{LF}','W_{LF} [kHz]','\mu_{LF}','E_{CS}','\sigma_{CS}','E_{S0}','L_{\epsilon}','\beta_{BB}'};
tLS = {'LOC','SOC'};
l = length(VarNames);
m = length(indq);
n = length(indLS);
rq1_md = zeros(1,m);
for k = [1,8]
        figure('Position',[95 90 616 591],'Color','w');
        for i = 1:m
            ax1 = subplot(m,1,i);
            
            for j = 1:n
            ind = indOK & indq{i} & indLS{j};   
            
            % Median values in a small interval
            for ii = 1:length(rc)
                indr = d.rhoc>rl(ii) & d.rhoc<rr(ii);
                indR = ind & indr;
                R_temp = d.(VarNames{k})(indR);
                if ~isempty(R_temp)
                    R(ii,i) = median(R_temp);
                    indN=(R_temp-median(R_temp))<0;
                    indP=(R_temp-median(R_temp))>0;
                    STD_median_N(ii,i) = abs(mean(R_temp(indN)-median(R_temp)));
                    STD_median_P(ii,i) = abs(mean(R_temp(indP)-median(R_temp)));
                else
                    R(ii,i) = NaN;
                    STD_median_N(ii,i) = NaN;
                    STD_median_P(ii,i) = NaN;
                end
            end
            
            % Median values within subintervals
            for ii = 1:length(rc)                
                VARdr = d.(VarNames{k})(ind & inddr{ii});
                if ~isempty(VARdr)
                    VAR.(VarNames{k})(ii,sub2ind([m,n],i,j)) = median(VARdr);
                else
                    VAR.(VarNames{k})(ii,sub2ind([m,n],i,j)) = NaN;
                end
            end 
            
            
           
%             ax1.Position(1)=0.05;
%             ax1.Position(3) = 0.43;
            %ax1.Position(4) = 0.30;
            hold on
            
            %ax1 = plot(d.rhoc(ind),d.(indF{k})(ind),'.','Color','c');
           % ax2 = plot((rc),R(:,i),'rs','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r');    
            
           % plot(d.rhoc(ind),d.(VarNames{k})(ind),'.','Color',[0 1 1]);
           if j==1
            ax2 = errorbar(rc,R(:,i),STD_median_N(:,i),STD_median_P(:,i),'^--',...
                'Color','r','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','r',...
                'MarkerFaceColor','r');
           elseif j==2
            ax3 = errorbar(rc,R(:,i),STD_median_N(:,i),STD_median_P(:,i),'s--',...
                'Color','b','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','b',...
                'MarkerFaceColor','b');
           end
            
            rq1 = mean(1./d.qpsi(ind));
            rq1_md(i) = rq1;
            if j==1
                 l1 = line([-rq1 -rq1],[0 0.4]*YLim{k}(2));
                 l1.Color = 'r'; lloc=l1;
            elseif j==2
                  l1 = line([-rq1 -rq1],[0.4 1]*YLim{k}(2));
                 l1.Color = 'b'; lsoc=l1;
            end
            l1.LineWidth = 2;
            l1.LineStyle = '-.';

            if j==1
                l2 = line([rq1 rq1],[0 0.4]*YLim{k}(2));
                l2.Color = 'r'; 
            elseif j==2
                l2 = line([rq1 rq1],[0.4 1]*YLim{k}(2));
                l2.Color = 'b'; 
            end
            l2.LineWidth = 2;
            l2.LineStyle = '-.';

            hold off
            
            ax = gca;
            ax.LineWidth=2;
            ax.Box = 'on';
            ax.XLim = [-1 0.6];
            ax.TickDir = 'out';
            ax.YLim = YLim{k};
            ax.FontSize = 16;
            ax.Position(4) = 0.19;
            
            if i==m
                ax.XTickLabel = {'-1','-0.5','0','0.5','1'};
                xlabel('\rho','FontSize',16);              
            else
                ax.XTickLabel = [];
            end
            
            %if i==round(m/2) && j==1;  ylabel(tVar{k},'FontSize',16); end
            
%            if i==1 
%                title(['Ohmic',', ',tVar{k},', Taylor'],'FontSize',12);
%            end
                     
           
            end
            
            t = text(-0,YLim{k}(2)*0.95,sprintf('%d < q_{\\psi} < %d',q(i),q(i+1)));
            t.FontSize = 14;
            t.BackgroundColor = 'white';
            t.VerticalAlignment = 'top';
            t.HorizontalAlignment = 'left';  
            
            t2=text(-0.95,YLim{k}(2)*0.95,txt{i},'FontSize',18);
        end
        lgd = legend([ax2 ax3 lloc lsoc],'LOC','SOC','q=1 for LOC','q=1 for SOC');
        lgd.Position = [0.4382    0.0073    0.1950    0.0463];
        lgd.Orientation = 'horizontal';
        lgd.FontSize = 14;

xlabel('\rho')
ylabel(tVar(k))

end