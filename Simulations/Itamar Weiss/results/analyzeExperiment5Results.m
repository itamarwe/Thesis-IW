close all;
clear all;

%Analyze the Experiment 1 Results

dpdUnknownSigRows = 1:10;
%dpdKnownSigRows = 10:18;
convUnknownSigRows = 11:20;
crbRows = 1:10;

snrCol = 11;
rmsXYCol = 19;
rmsVxVyCol = 20;
crbXYCol = 25;
crbVxVyCol = 26;


%Read the data from the file
[file dir] = uigetfile('*.xlsx');
[data,txt]= xlsread([dir file]);

%Plot the Position Performance Graph
figure;
semilogy(   data(dpdUnknownSigRows,snrCol),data(dpdUnknownSigRows,rmsXYCol),'b-s',...
            data(convUnknownSigRows,snrCol),data(convUnknownSigRows,rmsXYCol),'r-o');
hold on;
semilogy(data(crbRows,snrCol),data(crbRows,crbXYCol),'k-');
hold off;
legend('DPD Unknown Signals',...
    'Conventional Unknown Signals','CRB Known Signals');
title('Position Estimation Preformance - Circular - Random Signals','interpreter','latex','fontsize',16);
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:15);
set(gca,'xticklabel',-25:5:15);

%Plot the Velocity Performance Graph
figure;
semilogy(   data(dpdUnknownSigRows,snrCol),data(dpdUnknownSigRows,rmsVxVyCol),'b-s',...
            data(convUnknownSigRows,snrCol),data(convUnknownSigRows,rmsVxVyCol),'r-o');
hold on;
semilogy(data(crbRows,snrCol),data(crbRows,crbVxVyCol),'k-');
hold off;
legend('DPD Unknown Signals',...
    'Conventional Unknown Signals','CRB Known Signals');
title('Velocity Estimation Preformance - Circular Array - Random Signals','interpreter','latex','fontsize',16);
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:15);
set(gca,'xticklabel',-25:5:15);
