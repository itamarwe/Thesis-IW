function [tau,doppler]  = high_per_peak_est(CAF,Fs,tau_vec,dopp_vec);

% function [tau,doppler]  = high_per_peak_est(CAF,Fs,tau_vec,dopp_vec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate Doppler Shift and TOA  

% Inputs:
%-------------
% CAF                                       cross-ambiguity surface (Tau X Doppler)
% Fs                                           sampling rate
% tau_vec                                 vector of time-delays used in computing surface (1xTau)
% dopp_vec                             vector of Doppler frequenices used in computing surface (1xDoppler)

% Outputs:
%---------------
% tau                                         Estimated doppler shifts (step 1) 
% doppler                                 Estimated TOA  (step 1) 

% Author
%--------------
% Daniel Avergun (11.11.08)

%History
%--------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define delta values
delta_t = tau_vec(2) - tau_vec(1);
delta_f = dopp_vec(2) - dopp_vec(1);

% Get peak indices
[j1,dopp_pk_index] = max(max(CAF));
[j2,tau_pk_index] = max(max(CAF,[],2));

if dopp_pk_index ~= 1 && dopp_pk_index ~= size(CAF,2)
    % Get CAF values to fit
    dopp_values = CAF(tau_pk_index,dopp_pk_index-1:dopp_pk_index+1);
    % Compute the offset
    d_dopp = (dopp_values(1)-dopp_values(3))*delta_f/2/(dopp_values(1)+dopp_values(3)-2*dopp_values(2));
    % Calculate Tau and Doppler
    doppler = dopp_vec(dopp_pk_index) + d_dopp;
else
    doppler = dopp_vec(dopp_pk_index);
end

if tau_pk_index ~= 1 && tau_pk_index ~= size(CAF,1)
    % Get CAF values to fit
    tau_values = CAF(tau_pk_index-1:tau_pk_index+1,dopp_pk_index);
    % Compute the offset
    d_tau = (tau_values(1)-tau_values(3))*delta_t/2/(tau_values(1)+tau_values(3)-2*tau_values(2));
    % Calculate Tau and Doppler
    tau = tau_vec(tau_pk_index) + d_tau;
else
     tau = tau_vec(tau_pk_index);
end

