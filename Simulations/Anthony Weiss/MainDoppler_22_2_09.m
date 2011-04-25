clear all;
close all;
clc;
tic
CF1=0; %show cost function
CF11=0; %show cost function
DPDRun = 1;
STRun = 1;
Nexp = 50; % Number of Monte-Carlo trials 
Known_Unknown_Sig   = 0; % 1- Known Signal, 0 - Unknown Signals
c = 3e8;                                               % Speed of light - 3*10^2 m/s
Q = 1;                                                  % Number of emitters
K = 10; %5                                         % Number of interception intervals
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
NOF = 201;

%amplitudes = ones(K,NOF);
phases = rand(K,NOF)*0;
if NOF>1
    f0 = (1:NOF)./T0;
    f1 = f0(f0<0.4*Fs);
else
    f1 = 1/T0;
end

sss = ['V = ',int2str(ReceiverSpeed),' [m/s], f_c = ',num2str(fc/1e9),...
    ' [GHz], T_0 = ',int2str(round(T0/1e-6)),' [\mu Sec], F_s = ',num2str(Fs/1e6,'%3.1f'),' [MHz], BW = ',num2str(f1(end)/1e6,'%3.1f'), ' [MHz]'];
disp(sss)

%amplitudes = 1 + 0.3*rand(K,length(f1));
 amplitudes = sin(30*pi*f1/f1(end))./(30*pi*f1/f1(end));
b_vec=ones(1,K,L);
int = amplitudes(1,:)*cos( 2*pi*(f1')*(time_vec-T0/2)) ; 

figure
    plot(time_vec,int)
    drawnow
%keyboard

x_hat_DPD=zeros(length(SNRdB),Nexp);
y_hat_DPD=zeros(length(SNRdB),Nexp);
x_hat_ST=zeros(length(SNRdB),Nexp);
y_hat_ST=zeros(length(SNRdB),Nexp);
CRBp_Std = zeros(length(SNRdB),1);

[TOA,Doppler_shift] = Create_Tau_Doppler(EmitterPos,ReceiverPos,ReceiverVel,fc,c,L,K);
for snr_ind = 1:length(SNRdB)
    sigma2  = sigma2Vec(snr_ind);
    for nexp = 1:Nexp
        ttCPU = cputime;
        Delayed_sig = zeros(Ns,K,L);
        for k = 1:K
             for ell = 1:L
                DT = TOA(ell,k) - TOA(1,k);
                DD = Doppler_shift(ell,k) - Doppler_shift(1,k);
                int = amplitudes*cos( 2*pi*(f1')*(time_vec-DT-T0/2)) ; 
                %int=int./norm(int);
                Delayed_sig(:,k,ell)=(int.').*exp(1i*2*pi*DD*time_vec');
                E_sig = Delayed_sig(:,k,ell)'*Delayed_sig(:,k,ell)/Ns;
                sigma_sq = E_sig*sigma2;
                noise = sqrt(sigma_sq/2)*(randn(Ns,1)+1i*randn(Ns,1));
%                 if k==K 
%                     noise=noise*100;
%                 end
                r_vec_t(:,k,ell) =  Delayed_sig(:,k,ell) + noise ;  % (Eq.1 - new)
            end;%ell
        end;%k

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GRID SEARCH PARAMETERS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dLR_DPD =  2000;
        dLR_ST = dLR_DPD;
        dC_DPD = 1*dist_mult;
        dC_ST = dC_DPD;
        Xsearch_DPD = [-dLR_DPD,-dLR_DPD/2,0,dLR_DPD/2,dLR_DPD]+EmitterPos(1)+(rand-.5)*10;
        Ysearch_DPD = [-dLR_DPD,-dLR_DPD/2,0,dLR_DPD/2,dLR_DPD]+EmitterPos(2)+(rand-.5)*10;
        Xsearch_ST = Xsearch_DPD;
        Ysearch_ST = Xsearch_DPD;
        LRDiv_DPD = [2 10 1 ];
        CDiv_DPD = [10 10 1 ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DIRECT EMITTER LOCALIZATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if DPDRun
            %disp('DPD Search Started...')
            counter = 0;     
            while dLR_DPD > .1
                counter = counter+1;
                if CF1==1
                    Xsearch_DPD = (-150:10:150) +EmitterPos(1);
                    Ysearch_DPD = (-150:10:150)+EmitterPos(2);
                end
                
                [EmitPosEst_DPD,Xest_DPD,Yest_DPD,CF_DPD]  = DPDPosEst_f(time_vec,Xsearch_DPD,Ysearch_DPD,ReceiverPos,ReceiverVel,K,L,...
                                                           Q,fc,c,r_vec_t,Known_Unknown_Sig,CF1);

                dLR_DPD=dLR_DPD*.8;
                clear Xsearch_DPD Ysearch_DPD
                Xsearch_DPD = [-dLR_DPD,0,dLR_DPD]+Xest_DPD+(rand-.5)*0;
                Ysearch_DPD = [-dLR_DPD,0,dLR_DPD]+Yest_DPD+(rand-.5)*0;
            end

            x_hat_DPD(snr_ind,nexp)  = Xest_DPD;
            y_hat_DPD(snr_ind,nexp)  = Yest_DPD;

            Norm_mat(snr_ind,nexp) = norm(EmitPosEst_DPD - EmitterPos');
            disp(['SNR = ',num2str(SNRdB(snr_ind)),',  Exp_No = ', num2str(nexp)]);
                disp(['DPD_error = ',num2str( Norm_mat(snr_ind,nexp))]);
        end %DPD Run

        if STRun
            %[CRBp,CRB_f,CRB_tau] = crb(ReceiverPos,  ReceiverVel,  EmitterPos, fc,  f1, time_vec, amplitudes,  sigma_sq, TOA, Doppler_shift, K, c);
            tempo = ReceiverPos(:,:,1)-ReceiverPos(:,:,2);
            RANGE = sqrt(sum(tempo.*tempo,2));
            max_td = RANGE/c;
            [TD,FD] = TD_FD(r_vec_t,max_dop,max_td,time_vec);

            disp('Stein - Search Started...')
            if CF11==1
                    Xsearch_ST = (-150:10:150)*10 +EmitterPos(1);
                    Ysearch_ST = (-150:10:150)*10 +EmitterPos(2);
            end
            counter = 0;     
            while dLR_ST > .1
                counter = counter+1;
                
                [EmitPosEst_ST,Xest_ST,Yest_ST,CF_ST]  =  STPosEst(FD,TD,Xsearch_ST,...
                                                                    Ysearch_ST,K,L,ReceiverPos,ReceiverVel,fc,c,max(max_dop)^2,max(max_td)^2, CF11);

                dLR_ST = dLR_ST*.8;
                Xsearch_ST = Xest_ST + [-dLR_ST,0,dLR_ST];
                Ysearch_ST = Yest_ST + [-dLR_ST,0,dLR_ST];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The Estimated Postion
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x_hat_ST(snr_ind,nexp)  = Xest_ST;
            y_hat_ST(snr_ind,nexp)  = Yest_ST;

            Norm_mat_ST(snr_ind,nexp) = norm(EmitPosEst_ST - EmitterPos');
            
            disp(['ST Error = ',num2str( Norm_mat_ST(snr_ind,nexp))])
            disp(' ')
        end %STRun
    end% Nexp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CRAMER RAO LOWER BOUND
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     [CRBp,J] =   CRLB_d_9_2_09(Delayed_sig,amplitudes,phases,b_vec,sigma2,EmitterPos,ReceiverPos,ReceiverVel,Ni_vec,c,Nsc,K,L,Q,...
%         fc,Fs,max_dop,freq_vec,time_vec,Known_Unknown_Sig,TOA,Doppler_shift,f1);
[CRBp,CRB_f,CRB_tau] = crb(ReceiverPos,  ReceiverVel,  EmitterPos, fc,  f1, time_vec, amplitudes,  sigma_sq, TOA, Doppler_shift, K, c, T0);
     CRBp_Std(snr_ind)  =  sqrt(CRBp(1,1)+CRBp(2,2));

%         % For test only 12.2.09
        disp(' ')
         disp(['CRBp_Std = ',num2str(CRBp_Std(snr_ind))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate RMSE and Bias 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     %DPD
    if DPDRun  
        Bias_vec_DPD = [mean(x_hat_DPD(snr_ind,1:Nexp)); mean(y_hat_DPD(snr_ind,1:Nexp))] - EmitterPos';
        ErrorMat_DPD = [x_hat_DPD(snr_ind,1:Nexp);y_hat_DPD(snr_ind,1:Nexp)] - EmitterPos'*ones(1,Nexp);
        Pt_DPD = 1/Nexp * ErrorMat_DPD * ErrorMat_DPD';
        Poserr_DPD(1,snr_ind) = sqrt(Pt_DPD(1,1)+Pt_DPD(2,2));
        Bias_DPD(1,snr_ind) = norm([mean(ErrorMat_DPD(1,:)), mean(ErrorMat_DPD(2,:))]);

%             disp(['SNR,Exp ', '[' ,num2str(SNRdB(snr_ind)), ']' ,' Direct: Error ',num2str(Poserr_DPD(ns_ind,snr_ind)),...
%                 ' Direct: Bias ',num2str(Bias_DPD(ns_ind,snr_ind)), 'time:',num2str(cputime - ttCPU)]);
        disp(['SNR,Exp ', '[' ,num2str(SNRdB(snr_ind)), ']' ,' Direct: Error ',num2str(Poserr_DPD(1,snr_ind)),...
            ' Direct: Bias ',num2str(Bias_DPD(1,snr_ind)), 'CRBp_Std ', ...
             'time:',num2str(cputime - ttCPU)]);
    end %DPD 

    %Stein
    if STRun
        Bias_vec_ST = [mean(x_hat_ST(snr_ind,1:Nexp)); mean(y_hat_ST(snr_ind,1:Nexp))] - EmitterPos';
        ErrorMat_ST = [x_hat_ST(snr_ind,1:Nexp);y_hat_ST(snr_ind,1:Nexp)] - EmitterPos'*ones(1,Nexp);
        Pt_ST = 1/Nexp * ErrorMat_ST * ErrorMat_ST';
        Poserr_ST(1,snr_ind) = sqrt(Pt_ST(1,1) + Pt_ST(2,2));
        Bias_ST(1,snr_ind) = norm([mean(ErrorMat_ST(1,:)), mean(ErrorMat_ST(2,:))]);

%             disp(['SNR,Exp ', '[' ,num2str(SNRdB(snr_ind)), ']' ,' Stein: Error ',num2str(Poserr_ST(ns_ind,snr_ind)),...
%                 ' Stein: Bias ',num2str(Bias_ST(ns_ind,snr_ind)), 'time:',num2str(cputime - ttCPU)]);
         disp(['SNR,Exp ', '[' ,num2str(SNRdB(snr_ind)), ']' ,' Stein: Error ',num2str(Poserr_ST(1,snr_ind)),...
            ' Stein: Bias ',num2str(Bias_ST(1,snr_ind)), 'CRBp_Std ', ...
             'time:',num2str(cputime - ttCPU)]);
    end %Stein

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end;% SNR



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
if DPDRun && STRun
     semilogy(SNRdB,Poserr_DPD,'b-s','linewidth',1.2); hold on
     semilogy(SNRdB,Poserr_ST,'r-d','linewidth',1.2);
     semilogy(SNRdB,CRBp_Std,'k-','linewidth',1.2);hold off
     legend('DPD Algorithm','Stein Algorithm','CRB');
     title(sss)
elseif DPDRun
    semilogy(SNRdB,Poserr_DPD,'b-s','linewidth',1.2); hold on
    semilogy(SNRdB,CRBp_Std,'k-','linewidth',1.2);hold off
    legend('DPD Algorithm','CRLB');
elseif STRun
    semilogy(SNRdB,Poserr_ST,'b-s','linewidth',1.2); hold on
    semilogy(SNRdB,CRBp_Std,'k-','linewidth',1.2); hold off
    legend('Stein Algorithm','CRLB');
% else % For tests only  12.2.09
%     semilogy(SNRdB,CRBp_Std,'k-','linewidth',1.2)
%     legend('CRLB');
end

xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',SNRdB);
set(gca,'xticklabel',SNRdB);
toc
save 