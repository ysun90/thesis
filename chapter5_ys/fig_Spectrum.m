
close all

% Get spectrum
numc=319510;
[f,S]=spectrum_norm(numc);

% Plot the spectrum
figure('Pos',[388   147   698   526],'Color','w')
plot(f,10*log10(S),'b','LineWidth',2)
xlabel('Frequency [kHz]')
ylabel('Power spectrum [dB]');
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
ax.XTick=[-400 -200 0 200 400];

% Add text
text(10,-10,'\leftarrow CS','FontSize',20,'Color','r')
text(-150,-35,'LF \rightarrow','FontSize',20,'Color','r')
text(0,-52,'BB','FontSize',20,'Color','r','HorizontalAlignment','center')
text(400,-59,'N \downarrow','FontSize',20,'Color','r','HorizontalAlignment','right')
