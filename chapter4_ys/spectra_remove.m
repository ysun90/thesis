% Obtain spectra
[f1,S1]=spectrum_norm(2);
[f2,S2]=spectrum_norm(338011);
% plot spectra
figure('Color','w')
subplot(121)
plot(f1,lg(S1))
set(gca,'FontSize',18,'LineWidth',2)
xlabel('Frequency [kHz]')
ylabel('Power spectrum')
title('Spectrum with Low SNR')
text(-200,-30,'(a)','FontSize',20)

subplot(122)
plot(f2,lg(S2))
set(gca,'FontSize',18,'LineWidth',2)
xlabel('Frequency [kHz]')
ylabel('Power spectrum')
title('Spectrum with strong Doppler effect')
text(-200,-30,'(b)','FontSize',20)