%%
%close all
clear

%Load TS database
load('database_full.mat');
d = DatabaseYS;
d.Nla = d.Nl./(2*d.a); % [10^19 m-3]

% Load Taylor fit database
fit = load('FitTaylor.mat');
db = load('Database_v10_plasma_temp.mat');
d.rhoc = db.reflec.rhoc;
d.EBB = fit.EBB;
% cc=sqrt(d.Lepsi)./((3*10^8)./(d.f_plateau*10^9));
% d.EBB = fit.EBB./cc; % GHz
d.ELF = fit.ELF;
d.ECS = fit.ECS;
d.EN = fit.noise*1e3;
d.mu_BB = fit.mu_BB;
d.tL = fit.tL_BB;
d.K = fit.K_BB;
d.KtL = d.tL.*d.K;
d.WBB = fit.HWHM_BB;
d.Error = fit.resnorm;
clear DatabaseYS fit db

%Radil positions
dr = 0.05;
rc = dr-1:2*dr:1-dr;
rl = rc-dr;
ru = rc+dr;

% Filter database
indOK = d.rhoc>-1 & d.rhoc<1 & d.EquiTag==1 & d.Nla>0 & d.Error<0.05 & d.f_plateau<1e3 &... 
        ...%d.index_choc>0 &...             % New shots >40125
        ...%d.sigma_FB<250 & d.sigma_FB>30 &... % Remove noise signal
        ...%d.Resnorm<300 &...  % quantile(d.Resnorm,0.82) == Good fitting
        ...%d.alpha_FB<6 &...    % Remove noise signal or Doppler
        ...%d.SNR>5 &...   % Lower noise
        (d.EBB+d.ELF)>5*d.EN &...
        d.EBB<1 &...
        abs(d.mu_BB*180/pi)<50;   % Remove strong Doppler effect
% % Indexing for positive values
% indPOS = d.fV>0 & d.PowerFB>0 & d.alpha_FB>0 & d.IntF20>0 & d.IntS20>0 & d.IntBN>0;
% indOK = indOK & indPOS;
indOH = (d.PICRH + d.PECRH + d.PLH)<0.1;
% Indexing for qpsi
ql = [3,4,5,6];
qu = [4,5,6,10];
indq = cell(1,length(ql));
for i = 1:length(indq)
    indq{i} = d.qpsi>ql(i) & d.qpsi<qu(i);
end
% Indexing for factors
indF = {'EBB','WBB',};
% indF = {'EBB','tL','K','KtL','ELF','ECS','WBB'};

%% Plot turbulence properties vs radial position
tq = cell(1,length(ql));
for i = 1:length(tq)
    tq{i} = sprintf('%2.1f < q_{\\psi} < %2.1f',ql(i),qu(i));
end
YLim = {[0 1],[0 100]};
%YLim = {[0 1]/200,[0 100]};
YTitle = {'E_{BB}','w_{BB} [kHz]'};

l = length(indF);
m = length(indq);

rq1_md = zeros(1,m);
R = zeros(length(rc),m);
STD_median = zeros(length(rc),m);

for k = 1:1
        figure1=figure('Position',[95 90 616 591],'Color','w');
        for i = 1:m
            ind = indOK & indq{i} & indOH;
            
            % Median values in a small interval
            for ii = 1:length(rc)
                indr = d.rhoc>rl(ii) & d.rhoc<ru(ii);
                indR = ind & indr;
                R_temp = d.(indF{k})(indR);
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
            
            subplot(m,1,i);
            hold on
            ax1 = plot(d.rhoc(ind),d.(indF{k})(ind),'.','Color','c');
           % ax2 = plot((rc),R(:,i),'ko','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');    
            
            ax2 = errorbar(rc,R(:,i),STD_median_N(:,i),STD_median_P(:,i),'--o','Color','k','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
            %ax2 = plot(rc,R(:,i),'r+','MarkerSize',14);
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
            %ax.XLim = [-1.5 1.5];
            ax.YLim = YLim{k};
            ax.XLim=[-1 0.6];
            ax.TickDir = 'out';
            ax.FontSize = 20;
            ax.Position(4) = 0.17;
            ax.TickLength = [0.01 0.01];
            t = text(-0.97,YLim{k}(2)*0.95,sprintf('%d < q_{\\psi} < %d',ql(i),qu(i)));
            t.FontSize = 16;
            t.BackgroundColor = 'white';
            t.VerticalAlignment = 'top';
            t.HorizontalAlignment = 'left'; 
            if i==4
                ax.XTickLabel = {'-1.0','-0.5','0','0.5'};
                xlabel('\rho');
                ylabel(YTitle{k});
            else
                ax.XTickLabel = [];
            end  
        end
        %lgd = legend([ax2],'Median value');
        %lgd = legend([ax2],'Median value');
        %lgd = legend([ax1 ax2],'Original data','Median value');
        lgd = legend([ax1 ax2 l1],'Original data','Median value','q=1 surface');
        lgd.Position = [0.3672    0.9482    0.3006    0.0518];
        lgd.Orientation = 'horizontal';
        lgd.FontSize = 16;

end