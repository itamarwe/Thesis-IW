clear all; close all;

[data,txt]= xlsread('2010-06-24 Results Summary.xls','sheet1');

SNRCol = 12;
RMSXYCol = 20;
RMSVxVyCol = 21;
CRBXYCol = 25;
CRBVxVyCol = 26;

%**********************************8
%Performance Vs SNR
%**************************************
DPDUnknownWBRows = 1:6;
DPDKnownWBRows = 7:12;
DPDUnknownNBRows = 13:21;
DPDKnownNBRows = 22:30;

 
% Px:100 Py:300 Vx:200 Vy:200 Pulse Signal 101 samples Width  
%N = 512 Samples
%NB:Fs - 3000Hz B:1000Hz 
%WB:Fs:2^23 Hz B:0.4*2^23
%Circular Array With Radius 1000

figure;
semilogy(data(DPDUnknownWBRows,SNRCol),data(DPDUnknownWBRows,RMSXYCol),'b-s',...
    data(DPDKnownWBRows,SNRCol),data(DPDKnownWBRows,RMSXYCol),'r-o',...
data(DPDUnknownNBRows,SNRCol),data(DPDUnknownNBRows,RMSXYCol),'b-p',...
data(DPDKnownNBRows,SNRCol),data(DPDKnownNBRows,RMSXYCol),'r-v');
hold on;
semilogy(data(DPDKnownWBRows,SNRCol),data(DPDKnownWBRows,CRBXYCol),'k--');
semilogy(data(DPDKnownNBRows,SNRCol),data(DPDKnownNBRows,CRBXYCol),'k-');
hold off;
legend('DPD Unknown WB Signals','DPD known WB Signals','DPD Unknown NB Signals', 'DPD Known NB Signals',...
    'CRB Known WB Signals', 'CRB Known NB Signals');
title('Position Estimation Preformance - Circular Array','interpreter','latex','fontsize',16);
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:25);
set(gca,'xticklabel',-25:5:25);

figure;
semilogy(data(DPDUnknownWBRows,SNRCol),data(DPDUnknownWBRows,RMSVxVyCol),'b-s',...
    data(DPDKnownWBRows,SNRCol),data(DPDKnownWBRows,RMSVxVyCol),'r-o',...
data(DPDUnknownNBRows,SNRCol),data(DPDUnknownNBRows,RMSVxVyCol),'b-p',...
data(DPDKnownNBRows,SNRCol),data(DPDKnownNBRows,RMSVxVyCol),'r-v');
hold on;
semilogy(data(DPDKnownWBRows,SNRCol),data(DPDKnownWBRows,CRBVxVyCol),'k--');
semilogy(data(DPDKnownNBRows,SNRCol),data(DPDKnownNBRows,CRBVxVyCol),'k-');
hold off;
legend('DPD Unknown WB Signals','DPD known WB Signals','DPD Unknown NB Signals', 'DPD Known NB Signals',...
    'CRB Known WB Signals', 'CRB Known NB Signals');
title('Velocity Estimation Preformance - Circular Array','interpreter','latex','fontsize',16);
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m/s]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:25);
set(gca,'xticklabel',-25:5:25);



%**********************************8
%Performance Vs Circular Antena Elements
%**************************************
DPDUnknownWBRows = 31:38;
DPDKnownWBRows = 39:46;
%Antenna Elements Column
AECol = 2;
 


figure;
plot(data(DPDUnknownWBRows,AECol),data(DPDUnknownWBRows,RMSXYCol),'b-s',...
    data(DPDKnownWBRows,AECol),data(DPDKnownWBRows,RMSXYCol),'r-o');
hold on;
plot(data(DPDUnknownWBRows,AECol),data(DPDUnknownWBRows,CRBXYCol),'k-');
hold off;
legend('DPD Unknown WB Signals','DPD known WB Signals','CRB WB Known');
title('Position Estimation Preformance','interpreter','latex','fontsize',16);
xtext='Circular Array Antenna Elements';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',3:12);
set(gca,'xticklabel',3:12);

figure;
plot(data(DPDUnknownWBRows,AECol),data(DPDUnknownWBRows,RMSVxVyCol),'b-s',...
    data(DPDKnownWBRows,AECol),data(DPDKnownWBRows,RMSVxVyCol),'r-o');
hold on;
plot(data(DPDUnknownWBRows,AECol),data(DPDUnknownWBRows,CRBVxVyCol),'k-');
hold off;
legend('DPD Unknown WB Signals','DPD known WB Signals','CRB WB Known');
title('Velocity Estimation Preformance','interpreter','latex','fontsize',16);
xtext='Circular Array Antenna Elements';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m/s]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',3:12);
set(gca,'xticklabel',3:12);

%**********************************8
%Performance Vs Linear Antena Elements
%**************************************
DPDUnknownWBRows = 47:54;
DPDKnownWBRows = 55:62;
%Antenna Elements Column
AECol = 2;
 


figure;
plot(data(DPDUnknownWBRows,AECol),data(DPDUnknownWBRows,RMSXYCol),'b-s',...
    data(DPDKnownWBRows,AECol),data(DPDKnownWBRows,RMSXYCol),'r-o');
hold on;
plot(data(DPDUnknownWBRows,AECol),data(DPDUnknownWBRows,CRBXYCol),'k-');
hold off;
legend('DPD Unknown WB Signals','DPD known WB Signals','CRB Known WB');
title('Position Estimation Preformance','interpreter','latex','fontsize',16);
xtext='Linear Array Antenna Elements';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',3:12);
set(gca,'xticklabel',3:12);

figure;
plot(data(DPDUnknownWBRows,AECol),data(DPDUnknownWBRows,RMSVxVyCol),'b-s',...
    data(DPDKnownWBRows,AECol),data(DPDKnownWBRows,RMSVxVyCol),'r-o');
hold on;
plot(data(DPDKnownWBRows,AECol),data(DPDUnknownWBRows,CRBVxVyCol),'k-')
hold off;
legend('DPD Unknown WB Signals','DPD known WB Signals','CRB Known Signals');
title('Velocity Estimation Preformance','interpreter','latex','fontsize',16);
xtext='Linear Array Antenna Elements';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m/s]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',3:12);
set(gca,'xticklabel',3:12);


%**********************************8
%Performance Vs SNR - Linear Array
%**************************************
DPDUnknownWBRows = 63:71;
DPDKnownWBRows = 72:80;

 
% Px:100 Py:300 Vx:200 Vy:200 Pulse Signal 101 samples Width  
%N = 512 Samples
%WB:Fs:2^23 Hz B:0.4*2^23
%Linear Array With length 1000

figure;
semilogy(data(DPDUnknownWBRows,SNRCol),data(DPDUnknownWBRows,RMSXYCol),'b-s',...
    data(DPDKnownWBRows,SNRCol),data(DPDKnownWBRows,RMSXYCol),'r-o');
hold on;
semilogy(data(DPDKnownWBRows,SNRCol),data(DPDKnownWBRows,CRBXYCol),'k-');
hold off;
legend('DPD Unknown WB Signals','DPD known WB Signals','CRB Known Signals');
title('Position Estimation Preformance - Linear Array','interpreter','latex','fontsize',16);
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:25);
set(gca,'xticklabel',-25:5:25);

figure;
semilogy(data(DPDUnknownWBRows,SNRCol),data(DPDUnknownWBRows,RMSVxVyCol),'b-s',...
    data(DPDKnownWBRows,SNRCol),data(DPDKnownWBRows,RMSVxVyCol),'r-o');
hold on;
semilogy(data(DPDKnownWBRows,SNRCol),data(DPDUnknownWBRows,CRBVxVyCol),'k-');
hold off;
legend('DPD Unknown WB Signals','DPD known WB Signals','CRB Known Signals');
title('Velocity Estimation Preformance - Linear Array','interpreter','latex','fontsize',16);
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m/s]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:25);
set(gca,'xticklabel',-25:5:25);

