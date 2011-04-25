% Known Signals
clear all;
close all;
clc;
tic
CF1=0; %show cost function
CF11=0; %show cost function
DPDRun = 1;
STRun = 1;
Nexp = 10; % Number of Monte-Carlo trials 
Known_Unknown_Sig   = 0; % 1- Known Signal, 0 - Unknown Signals
c = 3e8;                                               % Speed of light - 3*10^2 m/s
Q = 1;                                                  % Number of emitters
K = 5; %5                                         % Number of interception intervals
L = 2;%3;                                              % Number of receivers
SNRdB =  (-15:10:-5) ;  %ajw
if CF1==1 || CF11==1
    SNRdB = 100;
end
ReceiverSpeed = 3e2;                      % receiver velocity (identical to all)
dist_mult = 1e3;
EmitterPos  = [4500 3700]; % [[1+5*rand,1+5*rand]*dist_mult;
fc  = 1e9;                             % carrier frequency
Ni_vec  =  zeros(1,K); %transmitter frequency shifts

max_dop = ReceiverSpeed/c*fc;
Fs=2^23;   % sampling low resolution
Ts = 1/Fs;
Ns = 2^9 ;   % - Number of samples in observation (after sampling the observation)
T0=Ns/Fs; % observation time


time_vec = (0:Ns-1)*Ts;
freq_vec = (0:Ns-1)*(Fs/Ns);
f_delta = freq_vec(2)-freq_vec(1);

ReceiverPos(:,1,1) =     (0:K-1)'*dist_mult ;
ReceiverPos(:,2,1) =    (0*ones(1,K))'*dist_mult;
ReceiverVel(:,:,1) =      (ones(K,1)*[ReceiverSpeed 0]);

ReceiverPos(:,1,2) =     (0:K-1)'*dist_mult +2000;
ReceiverPos(:,2,2) =     (0*ones(1,K))'*dist_mult +100;
ReceiverVel(:,:,2) =     (ones(K,1)*[ReceiverSpeed 0]);

figure
plot(ReceiverPos(:,1,1),ReceiverPos(:,2,1),'r->', ReceiverPos(:,1,2),ReceiverPos(:,2,2),'b->',...
    EmitterPos(1),EmitterPos(2),'rp') ; grid on
xlabel('Meters');
ylabel('Meters')
drawnow

%[TOA,Doppler_shift] = Create_Tau_Doppler(EmitterPos,ReceiverPos,ReceiverVel,fc,c,L,K);

sigma2Vec = 10.^(-SNRdB/10);                             % noise variance
sigma2dBVec = 10*log10(sigma2Vec);                              % noise variance [dB]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NOF = 201;

%amplitudes = ones(K,NOF);
%phases = rand(K,NOF)*0;
if NOF>1
    f0 = (1:NOF)./T0;
    f1 = f0(f0<0.4*Fs);
else
    f1 = 1/T0;
end


sss = ['V = ',int2str(ReceiverSpeed),' [m/s], f_c = ',num2str(fc/1e9),...
    ' [GHz], T_0 = ',int2str(round(T0/1e-6)),' [\mu Sec], F_s = ',...
    num2str(Fs/1e6,'%3.1f'),' [MHz], BW = ',num2str(f1(end)/1e6,'%3.1f'), ' [MHz], No. of Trials = ',int2str(Nexp)];
disp(sss)

%amplitudes = 1 + 0.3*rand(K,length(f1));
 amplitudes = sin(30*pi*f1/f1(end))./(30*pi*f1/f1(end));
b_vec=ones(1,K,L);
int = amplitudes(1,:)*cos( 2*pi*(f1')*(time_vec-T0/2)) ; 

figure
    plot(time_vec,int)
    xlabel('Time [\mu Sec]')
    ylabel('Amplitude')
    title('Signal Waveform')
    drawnow
%keyboard

[TOA,Doppler_shift] = Create_Tau_Doppler(EmitterPos,ReceiverPos,ReceiverVel,fc,c,L,K);

[Poserr_DPD_K,Poserr_ST_K,Bias_vec_DPD_K,Bias_vec_ST_K,CRBp_Std_K] = known_signals_routine(SNRdB,Nexp,Ns,K,L,Q,fc,c,T0,...
    time_vec,amplitudes,TOA,Doppler_shift,sigma2Vec,EmitterPos,ReceiverPos,ReceiverVel,CF1,CF11,f1,DPDRun,STRun,max_dop);

[Poserr_DPD_U,Poserr_ST_U,Bias_vec_DPD_U,Bias_vec_ST_U,CRBp_Std_U] = Unknown_signals_routine(SNRdB,Nexp,Ns,K,L,Q,fc,c,T0,time_vec,amplitudes,...
TOA,Doppler_shift,sigma2Vec,EmitterPos,ReceiverPos,ReceiverVel,CF1,CF11,f1,DPDRun,STRun,max_dop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File Name 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xtemp = clock;
ztemp = ceil(xtemp(end)*10);
if (Known_Unknown_Sig == 1)
    ResFile  =   '\KnSig';
else
    ResFile  =   '\UnknSig';
end;
curr_Dir = cd;
RunDate =   [curr_Dir,ResFile,num2str(ztemp),datestr(now,'ddmm'),'.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE RESULTS
%%%%%%%%%%%%%%%%

save (RunDate);
if DPDRun && STRun
%    save (RunDate,'b_vec','SNRdB','Poserr_DPD','Bias_DPD','Poserr_ST','Bias_ST','CRBp_Std','-append');
elseif DPDRun
%     save (RunDate,'SNRdB','Poserr_DPD','Bias_DPD','CRBp_Std','-append');
elseif STRun
%    save (RunDate,'b_vec','SNRdB','Poserr_ST','Bias_ST','CRBp_Std','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%
%PLOT Results
%%%%%%%%%%%%%%%%%%%

figure;
semilogy(SNRdB,Poserr_DPD_K,'b-s',SNRdB,Poserr_ST_K,'r-d',SNRdB,Poserr_DPD_U,'b-p',SNRdB,Poserr_ST_U,'r-h',SNRdB,CRBp_Std_K,'k-')
legend('DPD Known Sig','2-Step Known Sig','DPD Unknown Sig','2-Step Unknown Sig','CRB');
title(sss)
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',SNRdB);
set(gca,'xticklabel',SNRdB);
toc
save 