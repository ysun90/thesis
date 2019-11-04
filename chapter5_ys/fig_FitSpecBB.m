
% load data
load fig_FitSpecBB.mat

% reorganize data
f=S1.F;
s={S1.S,S2.S,S4.S};
BB={S1.BB,S2.BB,S4.BB};
LFs={S1.CS+S1.LF,S2.CS+S2.LF,S4.CS+S4.LF};
N={S1.N,S2.N,S4.N};
tt={'Gaussian','Exponential','Lorentzian'};

% plot
h=figure;
set(h,'Position',[129 177 1111 468],'Color','w')
for i=1:length(s)
   subplot(1,3,i)
   hold on
   p1=plot(f,lg(s{i}),'c-','LineWidth',4);
   p2=plot(f,lg(N{i}),'k-','LineWidth',2);
   p3=plot(f,lg(LFs{i}),'b-','LineWidth',2);
   f1=fill(f,lg(BB{i}),'g','FaceAlpha',0.5);
   p4=plot(f,lg(BB{i}),'r-','LineWidth',3);
   hold off
   xlabel('Frequency [kHz]')
   ylabel('Power spectrum [dB]')
   title(tt{i})
   legend([p1,p2,p3,p4],'Spectrum','N',...
       'LF + CS','BB')
   set(gca,'Box','on','FontSize',20,'LineWidth',2,...
       'XLim',[-450 450],'YLim',[-75 0])
   %text(0,-50,'E_{BB}','FontSize',20)
   annotation('textarrow',[.3 .5],[.3 .5],'String','E_{BB}',...
       'FontSize',16,'LineWidth',2)
end


