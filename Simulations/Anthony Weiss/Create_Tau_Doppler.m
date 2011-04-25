function [TOA,Doppler_shift] = Create_Tau_Doppler(EmitterPos,ReceiverPos,ReceiverVel,fc,c,L,K)

% function [TOA,Doppler_shift] = Create_Tau_Doppler(EmitterPos,ReceiverPos,ReceiverVel,fc,c,L,K,Q,Ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%-----------
% EmitterPos                      The positions of the emitters
% ReceiverPos                    The positions of the receivers
% ReceiverVel                      The velocities of the receivers
% K                                        Number of interception intervals
% L                                         Number of receivers
%Q                                         Number of emitters
%fc                                          carrier frequency
% c                                           Speed of light - 3*10^2 m/s
%Ns                                        Number of samples
 
% Outputs:
%--------------
% Doppler_shift                   Calculated doppler shifts (L+1,K)         
% TOA                                    Calculated TOA's (L+1,K)

%Author
%------------
% Daniel Avergun (6.12.08)
%
%History
%-----------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance = zeros(L,K);
TOA = zeros(L,K);
% DelaySam = zeros(L,K);
EmitterRelPos = zeros(2,L,K);
Vrel = zeros(L,K);
Doppler_shift = zeros(L,K);
for k = 1:K
    for ell = 1:L
        distance(ell,k) = norm(   ReceiverPos(k,:,ell)  -   EmitterPos  );
        TOA(ell,k) = distance(ell,k) /c;
        %     DelaySam(ell,k) = round(TOA(ell,k) /Ts);
        EmitterRelPos(:,ell,k) = [EmitterPos(1,1) - ReceiverPos(k,1,ell), EmitterPos(1,2) - ReceiverPos(k,2,ell)];
        Vrel(ell,k) = ReceiverVel(k,:,ell)*(EmitterRelPos(:,ell,k)./norm(EmitterRelPos(:,ell,k),2));
        Doppler_shift(ell,k) = fc*Vrel(ell,k)/c; %  (Eq.3- new)
    end%ell;
end%k;
TOA(3,:) = diff(TOA);
Doppler_shift(3,:) = diff(Doppler_shift);