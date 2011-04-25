function [CRB,CRB_f,CRB_tau] = crb(ReceiverPos,  ReceiverVel,  EmitterPos, fc,  f1, time_vec, amplitudes,  sigma_sq, TOA, Doppler_shift, K, c, T0)
CRB_f=zeros(K,1);
CRB_tau=CRB_f;
Ts=time_vec(2)-time_vec(1);
Ns=length(time_vec);
dm_dx = zeros(Ns,K);
dm_dy = dm_dx;
%dm_df1 = dm_dx;
%dm_dtau1 = dm_dx;
for k = 1:K
   DT = TOA(2,k) - TOA(1,k);
   DD = Doppler_shift(2,k) - Doppler_shift(1,k);
   int = amplitudes*cos( 2*pi*(f1')*(time_vec-DT-T0/2)) ; 
   %E_sig =int*int'/Ns;
   int_dot = -amplitudes*2*pi*diag(f1)*sin( 2*pi*(f1')*(time_vec-DT-T0/2)) ;
   s_tilde = int.';
   Ak = diag(exp(1i*2*pi*DD*time_vec));
   A_dot = 1i*2*pi*Ts*diag(exp(1i*2*pi*DD*time_vec).*(0:Ns-1));
   dm_df = A_dot*s_tilde;
   dm_dtau = -Ak*int_dot.';
   
   V1 = ReceiverVel(k,:,1); % row vec
   V2 = ReceiverVel(k,:,2); % row vec
   P1 = ReceiverPos(k,:,1);
   P2 = ReceiverPos(k,:,2);
   d1 = norm(P1 - EmitterPos);
   d2 = norm(P2 - EmitterPos);
     
   df_dx2 = V2(1)/d2 - V2*(EmitterPos - P2).'*(EmitterPos(1) - P2(1))/d2^3; 
   df_dx1 = V1(1)/d1 - V1*(EmitterPos - P1).'*(EmitterPos(1) - P1(1))/d1^3;
   df_dx = (fc/c)*( df_dx2 - df_dx1 );
   
   df_dy2 = V2(2)/d2 - V2*(EmitterPos - P2).'*(EmitterPos(2) - P2(2))/d2^3; 
   df_dy1 = V1(2)/d1 - V1*(EmitterPos - P1).'*(EmitterPos(2) - P1(2))/d1^3;
   df_dy = (fc/c)*( df_dy2 - df_dy1 );
   
   dtau_dx = (1/c)*( (EmitterPos(1) - P2(1))/d2 - (EmitterPos(1) - P1(1))/d1 ); 
   dtau_dy = (1/c)*( (EmitterPos(2) - P2(2))/d2 - (EmitterPos(2) - P1(2))/d1 );
   
   dm_dx(:,k) = dm_df*df_dx + dm_dtau*dtau_dx;
   dm_dy(:,k) = dm_df*df_dy + dm_dtau*dtau_dy;
   
   %dm_df1(:,k) = dm_df;
   %dm_dtau1(:,k) = dm_dtau;
   CRB_f(k) = sigma_sq/(2*real(dm_df'*dm_df));
   CRB_tau(k) = sigma_sq/(2*real(dm_dtau'*dm_dtau));
end
J(1,1) = 2*real( dm_dx(:)'*dm_dx(:));
J(2,2) = 2*real( dm_dy(:)'*dm_dy(:));
J(1,2) = 2*real( dm_dx(:)'*dm_dy(:));
J(2,1) = J(1,2);
CRB = sigma_sq*inv(J)*norm(int)^2;

   
   
   