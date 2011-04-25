function [Poserr_DPD,Poserr_ST,Bias_vec_DPD,Bias_vec_ST,CRBp_Std] = known_signals_routine(SNRdB,Nexp,Ns,K,L,Q,fc,c,T0,time_vec,amplitudes,...
TOA,Doppler_shift,sigma2Vec,EmitterPos,ReceiverPos,ReceiverVel,CF1,CF11,f1,DPDRun,STRun,max_dop,Delayed_sig,Delayed_sig0,sig)
x_hat_DPD=zeros(length(SNRdB),Nexp);
y_hat_DPD=zeros(length(SNRdB),Nexp);
x_hat_ST=zeros(length(SNRdB),Nexp);
y_hat_ST=zeros(length(SNRdB),Nexp);
CRBp_Std = zeros(length(SNRdB),1);
Norm_mat = zeros(length(SNRdB),Nexp);
Norm_mat_ST= Norm_mat;
Poserr_DPD = zeros(1,length(SNRdB));
Bias_DPD = Poserr_DPD;    Poserr_ST = Poserr_DPD ;   Bias_ST = Bias_DPD;
r_vec_t = zeros(Ns,K,L);
r_vec_t0 = r_vec_t;

disp('Generating Data (Known Signals)...')

for snr_ind = 1:length(SNRdB)
    sigma2  = sigma2Vec(snr_ind);
    sigma_sq = sigma2/Ns;
    [CRBp,CRB_f,CRB_tau] = crb_known_signals(ReceiverPos,  ReceiverVel,  EmitterPos, fc,  f1, time_vec, amplitudes,  sigma_sq, TOA, Doppler_shift, K, L, c, T0);
    CRB_f = CRB_f*(0.5+rand);
    CRB_tau = CRB_tau*(0.5+rand);
    for nexp = 1:Nexp
        ttCPU = cputime;
        for k = 1:K
            %sig(k,:) = amplitudes*cos( 2*pi*(f1')*(time_vec-T0/2)) ;
             for ell = 1:L
                noise = sqrt(sigma_sq/2)*(randn(Ns,1)+1i*randn(Ns,1));
                r_vec_t(:,k,ell) =  Delayed_sig(:,k,ell) + noise ;  % (Eq.1 - new)
                r_vec_t0(:,k,ell) =  Delayed_sig0(:,k,ell) + noise ;
            end;%ell
        end;%k
        
        dLR_DPD =  2000;
        dLR_ST = dLR_DPD;
     
        Xsearch_DPD = [-dLR_DPD,-dLR_DPD/2,0,dLR_DPD/2,dLR_DPD]+EmitterPos(1)+(rand-.5)*10;
        Ysearch_DPD = [-dLR_DPD,-dLR_DPD/2,0,dLR_DPD/2,dLR_DPD]+EmitterPos(2)+(rand-.5)*10;
        Xsearch_ST = Xsearch_DPD;
        Ysearch_ST = Xsearch_DPD;
   
        if DPDRun
            while dLR_DPD > .1
                if CF1==1
                    Xsearch_DPD = (-150:10:150) +EmitterPos(1);
                    Ysearch_DPD = (-150:10:150)+EmitterPos(2);
                end
                
%                 [EmitPosEst_DPD,Xest_DPD,Yest_DPD,CF_DPD]  = DPDPosEst_f(time_vec,Xsearch_DPD,Ysearch_DPD,...
%                    ReceiverPos,ReceiverVel,K,L,Q,fc,c,r_vec_t,Known_Unknown_Sig,CF1);
                [EmitPosEst_DPD,Xest_DPD,Yest_DPD]  = DPD_Known_Signals(time_vec,Xsearch_DPD,Ysearch_DPD,ReceiverPos,ReceiverVel,K,L,Q,fc,c,r_vec_t0,CF1,sig);

                dLR_DPD=dLR_DPD*.8;
                clear Xsearch_DPD Ysearch_DPD
                Xsearch_DPD = [-dLR_DPD,0,dLR_DPD]+Xest_DPD+(rand-.5)*0;
                Ysearch_DPD = [-dLR_DPD,0,dLR_DPD]+Yest_DPD+(rand-.5)*0;
            end

            x_hat_DPD(snr_ind,nexp)  = Xest_DPD;
            y_hat_DPD(snr_ind,nexp)  = Yest_DPD;

            Norm_mat(snr_ind,nexp) = norm(EmitPosEst_DPD - EmitterPos');
            disp('Known Signals')
            disp(['SNR = ',num2str(SNRdB(snr_ind)),',  Exp_No = ', num2str(nexp)]);
                disp(['DPD_error = ',num2str( Norm_mat(snr_ind,nexp))]);
        end %DPD Run

        if STRun
            %tempo = ReceiverPos(:,:,1)-ReceiverPos(:,:,2);
            %RANGE = sqrt(sum(tempo.*tempo,2));
            %max_td = RANGE/c;
            %[TD,FD] = TD_FD(r_vec_t,max_dop,max_td,time_vec);
            [TD,FD] = TD_FD_Known_Signals(r_vec_t0,max_dop,time_vec,sig);
%              TD,FD,
%              TOA(1:2,:).',Doppler_shift(1:2,:).'
%             keyboard
            disp('Stein - Search Started...')
            if CF11==1
                    Xsearch_ST = (-150:10:150)*10 +EmitterPos(1);
                    Ysearch_ST = (-150:10:150)*10 +EmitterPos(2);
            end
                 
            while dLR_ST > .1
                                
                [EmitPosEst_ST,Xest_ST,Yest_ST]  =  STPosEst0(FD,TD,Xsearch_ST,...
                          Ysearch_ST,K,L,ReceiverPos,ReceiverVel,fc,c,CF11,max(max(CRB_tau)),max(max(CRB_f)));

                dLR_ST = dLR_ST*.8;
                Xsearch_ST = Xest_ST + [-dLR_ST,0,dLR_ST];
                Ysearch_ST = Yest_ST + [-dLR_ST,0,dLR_ST];
            end

            x_hat_ST(snr_ind,nexp)  = Xest_ST;
            y_hat_ST(snr_ind,nexp)  = Yest_ST;

            Norm_mat_ST(snr_ind,nexp) = norm(EmitPosEst_ST - EmitterPos');
            
            disp(['ST Error = ',num2str( Norm_mat_ST(snr_ind,nexp))])
            disp(' ')
        end %STRun
    end% Nexp

 
    
    CRBp_Std(snr_ind)  =  sqrt(CRBp(1,1)+CRBp(2,2));

    disp(' ')
    disp(['CRBp_Std = ',num2str(CRBp_Std(snr_ind))]);
  
    if DPDRun  
        Bias_vec_DPD = [mean(x_hat_DPD(snr_ind,1:Nexp)); mean(y_hat_DPD(snr_ind,1:Nexp))] - EmitterPos';
        ErrorMat_DPD = [x_hat_DPD(snr_ind,1:Nexp);y_hat_DPD(snr_ind,1:Nexp)] - EmitterPos'*ones(1,Nexp);
        Pt_DPD = 1/Nexp * (ErrorMat_DPD * ErrorMat_DPD');
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
        Pt_ST = 1/Nexp * (ErrorMat_ST * ErrorMat_ST');
        Poserr_ST(1,snr_ind) = sqrt(Pt_ST(1,1) + Pt_ST(2,2));
        Bias_ST(1,snr_ind) = norm([mean(ErrorMat_ST(1,:)), mean(ErrorMat_ST(2,:))]);

%             disp(['SNR,Exp ', '[' ,num2str(SNRdB(snr_ind)), ']' ,' Stein: Error ',num2str(Poserr_ST(ns_ind,snr_ind)),...
%                 ' Stein: Bias ',num2str(Bias_ST(ns_ind,snr_ind)), 'time:',num2str(cputime - ttCPU)]);
         disp(['SNR,Exp ', '[' ,num2str(SNRdB(snr_ind)), ']' ,' Stein: Error ',num2str(Poserr_ST(1,snr_ind)),...
            ' Stein: Bias ',num2str(Bias_ST(1,snr_ind)), 'CRBp_Std ', ...
             'time:',num2str(cputime - ttCPU)]);
    end %Stein
end;% SNR