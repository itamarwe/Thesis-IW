close all;
clear all;

%Analyze the Experiment 1 Results

dpdUnknownSigRows = 1:9;
dpdKnownSigRows = 10:18;
convUnknownSigRows = 19:27;
crbRows = 1:9;

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
            data(dpdKnownSigRows,snrCol),data(dpdKnownSigRows,rmsXYCol),'r-o',...
            data(convUnknownSigRows,snrCol),data(convUnknownSigRows,rmsXYCol),'b-p');
hold on;
semilogy(data(crbRows,snrCol),data(crbRows,crbXYCol),'k-');
hold off;
legend('DPD Unknown Signals','DPD known Signals',...
    'Conventional Unknown Signals','CRB Known Signals');
title('Position Estimation Preformance - Linear Array','interpreter','latex','fontsize',16);
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:25);
set(gca,'xticklabel',-25:5:25);

%Plot the Velocity Performance Graph
figure;
semilogy(   data(dpdUnknownSigRows,snrCol),data(dpdUnknownSigRows,rmsVxVyCol),'b-s',...
            data(dpdKnownSigRows,snrCol),data(dpdKnownSigRows,rmsVxVyCol),'r-o',...
            data(convUnknownSigRows,snrCol),data(convUnknownSigRows,rmsVxVyCol),'b-p');
hold on;
semilogy(data(crbRows,snrCol),data(crbRows,crbVxVyCol),'k-');
hold off;
legend('DPD Unknown Signals','DPD known Signals',...
    'Conventional Unknown Signals','CRB Known Signals');
title('Velocity Estimation Preformance - Linear Array','interpreter','latex','fontsize',16);
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:25);
set(gca,'xticklabel',-25:5:25);
