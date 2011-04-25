function [CRB,CRB_f,CRB_tau] = crb_known_signals(ReceiverPos,  ReceiverVel,  EmitterPos, fc,  f1, time_vec, amplitudes,  sigma_sq, TOA, Doppler_shift, K, L, c, T0)
% Known Signals
CRB_f=zeros(K,L);
CRB_tau=CRB_f;
Ns=length(time_vec);
dm_dx = zeros(Ns,K);
dm_dy = dm_dx;
J=zeros(2,2);
for k = 1:K
    for ell = 1:L
       DD = Doppler_shift(ell,k);
       int = amplitudes*cos( 2*pi*(f1')*(time_vec-TOA(ell,k)-T0/2)) ; 
       int_dot = amplitudes*2*pi*diag(f1)*sin( 2*pi*(f1')*(time_vec-TOA(ell,k)-T0/2)) ;
       s_tilde = int.';
       A_ell_k = diag(exp(1i*2*pi*DD*time_vec));
       A_dot = 1i*2*pi*diag(exp(1i*2*pi*DD*time_vec).*time_vec);
       dm_df = A_dot*s_tilde;
       dm_dtau = A_ell_k*int_dot.';

       V1 = ReceiverVel(k,:,ell); % row vec
       P1 = ReceiverPos(k,:,ell);
       d1 = norm(P1 - EmitterPos);
       df_dx0 = V1(1)/d1 - V1*(EmitterPos - P1).'*(EmitterPos(1) - P1(1))/d1^3; 
       df_dx = (fc/c)*df_dx0;

       df_dy0 = V1(2)/d1 - V1*(EmitterPos - P1).'*(EmitterPos(2) - P1(2))/d1^3; 
       df_dy = (fc/c)*df_dy0;

       dtau_dx = (1/c)* (EmitterPos(1) - P1(1))/d1 ; 
       dtau_dy = (1/c)* (EmitterPos(2) - P1(2))/d1 ;

       dm_dx(:,k,ell) = dm_df*df_dx + dm_dtau*dtau_dx;
       dm_dy(:,k,ell) = dm_df*df_dy + dm_dtau*dtau_dy;

       CRB_f(k,ell) = sigma_sq/(2*real(dm_df'*dm_df));
       CRB_tau(k,ell) = sigma_sq/(2*real(dm_dtau'*dm_dtau));
       J(1,1) = J(1,1) + dm_dx(:,k,ell)'*dm_dx(:,k,ell);
       J(2,2) = J(2,2) + dm_dy(:,k,ell)'*dm_dy(:,k,ell);
       J(1,2) = J(1,2) + dm_dx(:,k,ell)'*dm_dy(:,k,ell);
    end
end
J(1,1) = 2*real( J(1,1));
J(2,2) = 2*real( J(2,2));
J(1,2) = 2*real( J(1,2));
J(2,1) = J(1,2);
CRB = sigma_sq*(J\eye(2))*norm(int)^2;

   
   
   