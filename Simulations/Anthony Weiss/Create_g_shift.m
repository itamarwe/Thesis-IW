function Doppler_shift = Create_g_shift(EmitterPos,ReceiverPos,ReceiverVel,fc,c,L,K,Q)

% function Doppler_shift = Create_g_shift(EmitterPos,ReceiverPos,ReceiverVel,f_delta,fc,c,L,K,Q,Ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%-----------
% EmitterPos                      The positions of the emitters
% ReceiverPos                    The positions of the receivers
% ReceiverVel                      The velocities of the receivers
% f_delta                              the size of a frequency bin
% K                                        Number of interception intervals
% L                                         Number of receivers
%Q                                         Number of emitters
%fc                                          carrier frequency
% c                                           Speed of light - 3*10^2 m/s
%Ns                                        Number of samples
 
% Outputs:
%--------------
% Doppler_shift                          

%Author
%------------
% Daniel Avergun (10.7.08)
%
%History
%-----------
% 1.9.08  Use round instead of ceil in the last calculation  (TW)
% 1.10.08 Return Doppler shift instead of g_shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EmitterRelPos = zeros(2,L,K,Q);
Vrel = zeros(L,K,Q);
Doppler_shift = zeros(L,K,Q);
g_shift = zeros(L,K,Q);

for ell = 1:L
    for k = 1:K
        for q = 1:Q
                EmitterRelPos(:,ell,k,q) = [EmitterPos(q,1) - ReceiverPos(k,1,ell), EmitterPos(q,2) - ReceiverPos(k,2,ell)];
                Vrel(ell,k,q) = ReceiverVel(k,:,ell)*(EmitterRelPos(:,ell,k,q)./norm(EmitterRelPos(:,ell,k,q),2)); 
                Doppler_shift(ell,k,q) = fc*Vrel(ell,k,q)/c; %  (Eq.3- new)             
%                 g_shift(ell,k,q) = round(Doppler_shift(ell,k,q)/f_delta);       
        end;%q
    end;%k
end;%ell