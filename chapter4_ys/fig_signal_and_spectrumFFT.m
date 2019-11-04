%[a,b,c]=ind2tab(320011)

% Load signal data
load #47189drefluc.mat
i = 11;
t = txc(:,i);
x = xc(:,i);

A = abs(x); % amplitude
phi = (angle(xc(:,i))); % phase
Acos = A.*cos(phi); % Acos(\phi)
Asin = A.*sin(phi); % Asin(\phi);

% Filter
%[b,a] = butter(6,0.99);
%x = filter(b,a,x);

t=(t-min(t))*1e3;
t(end)=10;


figure('Color','w');
subplot(4,2,1)
plot(t,A,'LineWidth',1.5);ylabel('Amplitude')
set(gca,'FontSize',18,'XTickLabel',[],'XLim',[min(t) max(t)]);
text(0,500,'(a)','FontSize',20)
subplot(4,2,3)  
plot(t,phi,'LineWidth',1.5);ylabel('Phase \phi')
set(gca,'FontSize',18,'XTickLabel',[],'XLim',[min(t) max(t)]);
text(0,5,'(b)','FontSize',20)
subplot(4,2,5)
plot(t,Acos,'LineWidth',1.5);ylabel('Acos\phi')
set(gca,'FontSize',18,'XTickLabel',[],'XLim',[min(t) max(t)]);
text(0,500,'(c)','FontSize',20)
subplot(4,2,7)
plot(t,Asin,'LineWidth',1.5);xlabel('t [ms]'); ylabel('Asin\phi');
set(gca,'FontSize',18,'XLim',[min(t) max(t)]);
text(0,500,'(d)','FontSize',20)

subplot(4,2,[2,4,6,8])
% FFT
nfft=1024;
[pxx,f] = pwelch(x,1024,512,nfft,1e6);
f = linspace(-max(f)/2,max(f)/2,nfft);
f = f'/1e3; % kHz
pxx = fftshift(pxx);
%figure('Color','w');
plot(f,lg(pxx),'LineWidth',2); 
set(gca,'FontSize',18,'LineWidth',2);
xlabel('Frequency [kHz]')
ylabel('Power spectrum [dB]')
title('FFT from A^{i\phi}')
text(0,20,'(e)','FontSize',20)