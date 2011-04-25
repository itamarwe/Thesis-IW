clear all;
close all;
clc;
tic
CF1=0; %show cost function
CF11=0; %show cost function
DPDRun = 1;
STRun = 1;
Nexp = 100; % Number of Monte-Carlo trials 
Known_Unknown_Sig   = 0; % 1- Known Signal, 0 - Unknown Signals
c = 3e8;                                               % Speed of light - 3*10^2 m/s
Q = 1;                                                  % Number of emitters
K = 5; %5                                         % Number of interception intervals
L = 2;%3;                                              % Number of receivers
SNRdB =  (-15:10:45) ;  %ajw
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

ReceiverPos(:,1,1) =     (0:K-1)'*dist_mult +1000;
ReceiverPos(:,2,1) =    (0*ones(1,K))'*dist_mult;
ReceiverVel(:,:,1) =      (ones(K,1)*[ReceiverSpeed 0]);

ReceiverPos(:,1,2) =     (0:K-1)'*dist_mult +3000;
ReceiverPos(:,2,2) =     (0*ones(1,K))'*dist_mult;
ReceiverVel(:,:,2) =     (ones(K,1)*[ReceiverSpeed 0]);

figure
plot(ReceiverPos(:,1,1),ReceiverPos(:,2,1),'r-o', ReceiverPos(:,1,2),ReceiverPos(:,2,2),'b-*',...
    EmitterPos(1),EmitterPos(2),'rp') ; grid on
drawnow

%[TOA,Doppler_shift] = Create_Tau_Doppler(EmitterPos,ReceiverPos,ReceiverVel,fc,c,L,K);

sigma2Vec = 10.^(-SNRdB/10);                             % noise variance
sigma2dBVec = 10*log10(sigma2Vec);                              % noise variance [dB]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Norm_mat = zeros(length(SNRdB),Nexp);
Norm_mat_ST= Norm_mat;
Poserr_DPD = zeros(1,length(SNRdB));
Bias_DPD = Poserr_DPD;    Poserr_ST = Poserr_DPD ;   Bias_ST = Bias_DPD;
r_vec_t = zeros(Ns,K,L);
r_vec_f = zeros(Ns,K,L);
r_vec_f_temp = zeros(2*Ns,K,L);

disp('Generating Data...')
NOF = 51;

%amplitudes = ones(K,NOF);
phases = rand(K,NOF)*0;
if NOF>1
    f0 = (1:NOF)./T0;
    f1 = f0(f0<0.4*Fs);
else
    f1 = 1/T0;
end

sss = ['V = ',int2str(ReceiverSpeed),' [m/s], f_c = ',num2str(fc/1e9),' [GHz], T_0 = ',int2str(round(T0/1e-6)),' [\mu Sec], F_s = ',num2str(Fs/1e6,'%3.1f'),' [MHz], BW = ',num2str(f1(end)/1e6,'%3.1f'), ' [MHz]'];
disp(sss)
