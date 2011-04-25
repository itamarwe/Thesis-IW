function [Est_Doppler_Shift,Est_TOA] =    Est_Doppler_toa2(r_vec_f_temp,r_vec_t,r_vec_f,Ns,Fs,Ts,K,time_vec,freq_vec,max_dop,max_TOA,Known_Unknown_Sig);

%function [Est_Doppler_Shift,Est_TOA] =
%Est_Doppler_toa2(r_vec_t,r_vec_f,Ns,Fs,Ts,K,time_vec,freq_vec,max_dop,max_TOA,Known_Unknown_Sig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate Doppler Shift and TOA  

% Inputs:
%-------------
% r_vec_f                                   the vector of the obsreve signal at the receivers in freq domain
% Ns                                            Number of samples
% Fs                                             sampling rate
% Ts                                             sampling interval
% K                                              Number of interception intervals
% freq_vec                                 frequencies vector
% time_vec                                 time vector
% max_dop                                maximum doppler shift ( = ReceiverSpeed/c*fc)
% max_TOA                              maximum TOA shift ( = max_distance/c)
%Known_Unknown_Sig        Know / Unknown signal flag

% Outputs:
%---------------
% Est_Doppler_Shift            Estimated doppler shifts (step 1) 
% Est_TOA                              Estimated TOA  (step 1) 

% Author
%--------------
% Daniel Avergun (23.10.08)

%History
%--------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Est_Doppler_Shift = zeros(1,K);
Est_TOA = zeros(1,K);

time_vec_temp = (0:2*Ns-1)*Ts;
freq_vec_temp = (0:(2*Ns-1))*(Fs/(2*Ns));
dopp_vec_temp = freq_vec_temp(2:end) - Fs/2;
    
T = Ts*Ns;            % total observation time (per k)
nof_time_bin = 1000;

% %configure grid parameters
TOA_vec = linspace(-max_TOA,max_TOA,nof_time_bin);
dopp_vec = freq_vec(2:end) - Fs/2;

for k = 1:K
    for toa_ind = 1:nof_time_bin
        FFT1_temp = r_vec_f_temp(:,k,1);
        FFT1 = FFT1_temp(1:Ns);
        FFT2_temp = r_vec_f_temp(:,k,2).*exp(+j*2*pi*TOA_vec(toa_ind)*freq_vec_temp');
        FFT2 = FFT2_temp(1:Ns);
        CAF_temp(toa_ind,:) =  abs(xcorr(FFT2,FFT1).' );
    end %toa

    % Get tau and Doppler measurements from surface peak
%     [tau_est_curv2,doppler_est_curv2] = curve_fit(CAF_temp3,TOA_vec,dopp_vec_temp);
    
    % High persission peak estimation
    
     [tau_est_high,doppler_est_high] = high_per_peak_est(CAF_temp,Fs,TOA_vec,dopp_vec_temp);
      
    clear('FFT1','FFT2','CAF_temp');
%     k
    Est_Doppler_Shift(k) = doppler_est_high;
    Est_TOA(k) = tau_est_high;
end %k

