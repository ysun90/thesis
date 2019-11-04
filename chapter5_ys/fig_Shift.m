load #44468CutOoff.mat

figure('Pos',[327    49   668   635])

subplot(311);
plot(equi.rhog(1,:).*equi.signeg(1,:),equi.ne_fit(1,:)/1e19,'b-','LineWidth',2); hold on
plot(equi.rhog(1,:).*equi.signeg(1,:),equi.ne_refluc_ref(1,:)/1e19,'r-.','LineWidth',2); 
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
xlim([-1 1])
ylim([0 3.5])
%xlabel('r/a');
ylabel('n_e [\times 10^{19} m^{-3}]');
legend('Interferometry','Reflectometry')
text(-0.75,3,'(a)','FontSize',18)

subplot(312)
plot(equi.rhog(1,:).*equi.signeg(1,:),equi.fxh(1,:)/1e9,'b-','LineWidth',2); hold on
plot(equi.rhog(1,:).*equi.signeg(1,:),equi.fxh_refluc_ref(1,:)/1e9,'r-.','LineWidth',2); 
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
xlim([-1 1])
ylim([100 150])
%xlabel('r/a');
ylabel('F_{xh} [GHz]');
legend('Interferometry','Reflectometry')
text(-0.75,120,'(b)','FontSize',18)

subplot(313);
drhoc=(test.Rc1_refluc_ref(1,:)-test.Rc1_interf(1,:))*100;
plot(test.rhoc1_interf(1,:),drhoc,'k.','MarkerSize',20); 
ax=gca;
ax.LineWidth=2;
ax.FontSize=20;
xlim([-1 1])
ylim([0 8])
xlabel('\rho ')
ylabel('Shift [cm]')
text(-0.75,6,'(c)','FontSize',18)