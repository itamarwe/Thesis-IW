function [EmitPosEst,Xest,Yest,CF] = DPDPosEst_f(time_vec,Xsearch,Ysearch,ReceiverPos,ReceiverVel,K,L,Q,fc,c,Nsc,r_vec_t,freq_vec,f_delta,Known_Unknown_Sig,CF1)

%function [EmitPosEst,Xest,Yest,CF] = DPDPosEst_f(time_vec,Xsearch,Ysearch,ReceiverPos,ReceiverVel,K,L,Q,fc,c,Ns,r_vec_t,freq_vec,f_delta,Known_Unknown_Sig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%-----------
% Xsearch, Ysearch               The coordinates of the grid search for the emitters
% ReceiverPos                        The positions of the receivers
% ReceiverVel                        The velocities of the receivers
% freq_vec                               frequencies vector
% time_vec                               time vector
% r_vec_t                                 the vector of the obsreved signal at the receivers 
% f_delta                                  the size of a frequency bin
% K                                            Number of interception intervals
% L                                             Number of receivers
%Q                                             Number of emitters
%fc                                             carrier frequency
% c                                              Speed of light - 3*10^2 m/s
%Ns                                            Number of samples
%Known_Unknown_Sig        Know / Unknown signal flag


% Outputs:
%--------------
% EmitPosEst                              the estimated position of the emitter 
% Xest, Yest                                  the coordinates of the estimated position
%CF                                                Cost function

%Author
%------------
% Daniel Avergun (10.7.08)
%
%History
%-----------
% 20.7.08   Known signal ability was added  % by DA
% 24.8.08   1) Downsample the signal in the Known Signal ability 
 %                2) Insert signal matrix in a time domain instead of freq domain
  % 1.9.08   1) In calculating F_mat the shift has to be left (-1) not right
 %               2) In calculating A_mat the sign in the exponent was fixed
 %                    to  -1     % by TW
 %11.9.08  More efficient run time changes were added 
 %                1) "F_mat" is not calculated - instead the circshift is performed on "Vec_temp"
 %                2) "A_mat"  is not calculated - instead "Vec_temp" is calculated via "A_vec"
 %1.10.08   New version of the calculation was applied (Taken from Tony's version from 15.9.08)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CF = zeros(length(Xsearch),length(Ysearch));
Tsc = time_vec(2)-time_vec(1);

for x_ind = 1:length(Xsearch)
       for y_ind = 1:length(Ysearch)
                EmitPostmp = [Xsearch(x_ind),Ysearch(y_ind)];
                Dop_shift = Create_g_shift(EmitPostmp,ReceiverPos,ReceiverVel,fc,c,L,K,Q);
                CF_temp = 0;
                for k = 1:K
%                     Vk = zeros(Ns/2,L);
                    %Vk = zeros(Nsc,L);
                    for ell = 1:L
                        TOA_rel = norm(ReceiverPos(k,:,ell) - EmitPostmp)/c; %TOA_rel = Relative distance /c
                        sig_t = r_vec_t(:,k,ell);
                        s1=floor(TOA_rel/Tsc);
                        s2=ceil(TOA_rel/Tsc);
                        T1=TOA_rel-s1*Tsc;
                        T2=Tsc-T1;
                        int=circshift(sig_t,-s1)*T2/Tsc + circshift(sig_t,-s2)*T1/Tsc;
                        %int=circshift(sig_t,-s1);
                        %sig_time = real(ifft(fftshift(fftshift(fft(sig_t)).*exp(i*2*pi*TOA_rel*(freq_vec-1/Ts/2)'  ) )   )); %time shift
                        ttt1 = int.*exp(-i*2*pi*Dop_shift(ell,k,1)*time_vec'); %freq shift
                        %ttt = downsample(ttt1,10);
                        FFT1 = fft(ttt1);   % Insert zero padding  length Ns
                        FFT_new = FFT1;
                        tt2 = FFT_new;
                        Vk(:,ell) = tt2;
            
                    end;%ell
                    if Known_Unknown_Sig % Known Signals
                        %
                    else % Unknown Signals
                        Lambda  = svd(Vk'*Vk); % (Eq.22 - old)
                        CF_temp =  CF_temp + abs(Lambda(1));
                    end %if
                end;%k
                CF(x_ind,y_ind) = CF_temp;
       end; % y_ind
end; % x_ind
if CF1
 figure
 mesh(CF), 
 keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the optimum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u1,v1] = find(CF == max(max(CF))); 
Xest = Xsearch(u1); 
Yest = Ysearch(v1);
EmitPosEst = [Xest,Yest]';
% figure,surf(CF), shading interp

