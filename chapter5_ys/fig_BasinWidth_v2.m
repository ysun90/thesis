%% Plot the half width of energy basin with r_q1
rq1_md=[0.2724    0.2174    0.1852    0.1319];
dr=0.01;
rBasinLeft = [-0.33 -0.29 -0.29 -0.23];
rBasinRight = [0.19 0.13 0.07 0.03];

figure('Color','w','Pos',[138         147        1111         496])

subplot(121)
hold on
errorbar(rq1_md,(rBasinRight-rBasinLeft)/2,ones(1,length(rq1_md))*dr*2,...
  'bo','MarkerSize',20);
plot([0.1 0.3],[0.1 0.3],'r-','LineWidth',3)
hold off
set(gca,'LineWidth',2,'box','on')
axis square
set(gca,'FontSize',20)
xlabel('\rho_{q=1}')
ylabel('Half width of the E_{BB} basin (w/2)')
text(0.1,0.3,'(a)','FontSize',20)

subplot(122)
hold on
errorbar(rq1_md,abs(rBasinLeft),ones(1,length(rq1_md))*dr*2,...
  'b<','MarkerSize',20);
errorbar(rq1_md,abs(rBasinRight),ones(1,length(rq1_md))*dr*2,...
  'b>','MarkerSize',20);
plot([0 0.4],[0 0.4],'r-','LineWidth',3)
hold off
set(gca,'LineWidth',2,'box','on')
axis square
set(gca,'FontSize',20)
xlabel('\rho_{q=1}')
ylabel('Width of the E_{BB} basin')
legend('HFS','LFS')
text(0.1,0.3,'(b)','FontSize',20)