function [CRBp,J] =   CRLB_d_9_2_09(Delayed_sig,amplitudes,phases,b_vec,sigma2,EmitterPos,ReceiverPos,ReceiverVel,Ni_vec,c,Ns,K,L,Q,...
                                                        fc,Fs,max_dop,freq_vec,time_vec,Known_Unknown_Sig,TOA,Doppler_shift,f1);

% function [CRBp,J] =   CRLB_d_9_2_09(Delayed_sig,amplitudes,phases,b_vec,sigma2,EmitterPos,ReceiverPos,ReceiverVel,Ni_vec,c,Ns,K,L,Q,...
%                                                         fc,Fs,max_dop,freq_vec,time_vec,Known_Unknown_Sig,TOA,Doppler_shift,f1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%-----------
% Delayed_sig
% amplitudes
% phases
% b_vec
% EmitterPos                          The positions of the emitters
% ReceiverPos                        The positions of the receivers
% ReceiverVel                         The velocities of the receivers
% Ni_vec                                   Transmited frequency shifts
% c                                              Speed of light - 3*10^2 m/s
%Ns                                           Number of samples
% K                                             Number of interception intervals
% L                                             Number of receivers
%Q                                             Number of emitters
%fc                                             carrier frequency
% Fs                                            sampling rate
% max_dop                              maximum doppler shift ( = ReceiverSpeed/c*fc)
% freq_vec                               frequencies vector
% Known_Unknown_Sig     Know / Unknown signal flag
% TOA
%Doppler_shift
% f1
 
% Outputs:
%--------------
% CRBp                                    Cramer Rao Boundry for the position estimation
% J                                             Fisher information matrix

%Author
%------------
% Daniel Avergun (1.2.09)
%
%History
%-----------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = length(freq_vec);
FFT_length = 2*Ns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calaculating Jpp (2Q x 2Q matrix)  (Q=1 in this case)   - 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jxx = zeros(Q,Q);
Jyy = zeros(Q,Q);
Jxy = zeros(Q,Q);
Jyx = zeros(Q,Q);

% Calculate  dx , dy and df derivatives of the signal 
%-------------------------------------------------------------------------

% calculate new TOA and Doppler shifts according to dx and dy
dx = 1e-3;
EmitterPosTemp_dX  = EmitterPos + ones(Q,1)*[dx,0];
[TOA_dX,Doppler_shift_dX] = Create_Tau_Doppler(EmitterPosTemp_dX,ReceiverPos,ReceiverVel,fc,c,L,K);

dy =  1e-3;
EmitterPosTemp_dY  = EmitterPos+ones(Q,1)*[0,dy];
[TOA_dY,Doppler_shift_dY] = Create_Tau_Doppler(EmitterPosTemp_dY,ReceiverPos,ReceiverVel,fc,c,L,K);

% set df
df = max_dop/20;

for ell =1:L
    for k = 1:K
        % Calculate the delayed signals respectively to the new TOA,
        % doppler and frequency shifts 
                    Delayed_sig_dX(:,k,ell) = amplitudes(k,:)*sin( 2*pi*(f1' + Doppler_shift_dX(ell,k) + Ni_vec(k))*...
                                                (time_vec - TOA_dX(ell,k)) + phases(k,:)'*ones(size(time_vec)));
                    Delayed_sig_dY(:,k,ell) = amplitudes(k,:)*sin( 2*pi*(f1' + Doppler_shift_dY(ell,k) + Ni_vec(k))*...
                                                (time_vec - TOA_dY(ell,k)) + phases(k,:)'*ones(size(time_vec)));
                    Delayed_sig_dF(:,k,ell) = amplitudes(k,:)*sin( 2*pi*(f1' + Doppler_shift(ell,k) + Ni_vec(k) + df)*...
                                                (time_vec - TOA(ell,k)) + phases(k,:)'*ones(size(time_vec)));
                                            
                    Delayed_sig_f_temp = fft(Delayed_sig(:,k,ell),FFT_length);
                    Delayed_sig_f(:,k,ell) =   Delayed_sig_f_temp(1:FFT_length/2);
                    Delayed_sig_dX_f_temp = fft(Delayed_sig_dX(:,k,ell),FFT_length);
                    Delayed_sig_dX_f(:,k,ell) =   Delayed_sig_dX_f_temp(1:FFT_length/2);
                    Delayed_sig_dY_f_temp = fft(Delayed_sig_dY(:,k,ell),FFT_length);
                    Delayed_sig_dY_f(:,k,ell) =   Delayed_sig_dY_f_temp(1:FFT_length/2);
                    Delayed_sig_dF_f_temp = fft(Delayed_sig_dF(:,k,ell),FFT_length);
                    Delayed_sig_dF_f(:,k,ell) =   Delayed_sig_dF_f_temp(1:FFT_length/2);
                    
                    % Calculate the derivatives 
                    temp_Diff_sig_dX(:,k,ell) =  [Delayed_sig_dX_f(:,k,ell) -  Delayed_sig_f(:,k,ell)]/dx;
                    temp_Diff_sig_dY(:,k,ell) =  [Delayed_sig_dY_f(:,k,ell) -  Delayed_sig_f(:,k,ell)]/dy;
                    temp_Diff_sig_dF(:,k,ell) =  [Delayed_sig_dF_f(:,k,ell) -  Delayed_sig_f(:,k,ell)]/df;
    end %k
end %ell
                                            

for ell = 1:L
    for k = 1:K

        temp1 = abs(b_vec(:,k,ell)).^2 * temp_Diff_sig_dX(:,k,ell)' *temp_Diff_sig_dX(:,k,ell);    
        temp2 = abs(b_vec(:,k,ell)).^2 * temp_Diff_sig_dY(:,k,ell)' *temp_Diff_sig_dY(:,k,ell);
        temp3 = abs(b_vec(:,k,ell))^2 * temp_Diff_sig_dX(:,k,ell)' *temp_Diff_sig_dY(:,k,ell);
                   
        Jxx = Jxx + 2/sigma2*real(temp1);       
        Jyy = Jyy + 2/sigma2*real(temp2);       
        Jxy = Jxy + 2/sigma2*real(temp3);        
        Jyx = Jxy;        
    end;%k
end;%ell

Jpp =   [Jxx Jxy;Jyx Jyy];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calaculating Jpf (2Q x K matrix)       - 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jxf = zeros(Q,K);
Jyf = zeros(Q,K);

for k   =   1:K
    for ell =   1:L
        temp1 = abs(b_vec(:,k,ell))^2 * temp_Diff_sig_dX(:,k,ell)' *temp_Diff_sig_dF(:,k,ell);
        temp2 = abs(b_vec(:,k,ell))^2 * temp_Diff_sig_dY(:,k,ell)' *temp_Diff_sig_dF(:,k,ell);
        Jxf(k) = Jxf(k) + 2/sigma2*real(temp1);       
        Jyf(k) = Jyf(k) + 2/sigma2*real(temp2);       
    end;%ell
end;%k

Jfx = Jxf;
Jfy = Jyf;

Jpf = [Jxf;Jyf];
Jfp = Jpf';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calaculating Jpbar(b) (2Q x KL matrix)  & Jptilde(b) (2Q x KL matrix)   - 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jxbarb = zeros(Q,K*L);
Jxtildeb = zeros(Q,K*L);
Jybarb = zeros(Q,K*L);
Jytildeb = zeros(Q,K*L);

for k = 1:K
    for ell = 1:L
        temp1 = conj(b_vec(:,k,ell)) * temp_Diff_sig_dX(:,k,ell)' *Delayed_sig_f(:,k,ell);
        temp2 = conj(b_vec(:,k,ell)) * temp_Diff_sig_dY(:,k,ell)' *Delayed_sig_f(:,k,ell);
        Jxbarb(1,k+K*(ell-1)) = 2/sigma2*real(temp1);       
        Jybarb(1,k+K*(ell-1)) = 2/sigma2*real(temp2);
        Jxtildeb(1,k+K*(ell-1) ) = -2/sigma2*imag(temp1);
        Jytildeb(1,k+K*(ell-1) ) = -2/sigma2*imag(temp2);
    end;%ell
end;%k

Jbarbx =  Jxbarb';
Jtildebx  = Jxtildeb';
Jbarby = Jybarb';
Jtildeby  = Jytildeb';

Jpb = [Jxbarb Jxtildeb;Jybarb Jytildeb];
Jbp = [Jbarbx Jbarby;Jtildebx Jtildeby];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calaculating Jff (K x K matrix)  &  - 2
%Jfbar(b) (K x KL matrix) & Jftilde(b)  (K x KL matrix)   -   6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jff = zeros(K,K);
for k = 1:K 
    temp  =  0;
    for ell = 1:L
        temp = temp +  2/sigma2*real(abs(b_vec(:,k,ell))^2 * temp_Diff_sig_dF(:,k,ell)' *temp_Diff_sig_dF(:,k,ell));       
    end;%ell
    Jff(k,k)  = temp;
end;%k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jfbarb = zeros(K,K*L);
Jftildeb = zeros(K,K*L);
for k  = 1:K
    for ell = 1:L
                temp = conj(b_vec(:,k,ell))*temp_Diff_sig_dF(:,k,ell)' *Delayed_sig_f(:,k,ell);
                Jfbarb(k,k+K*(ell-1) )  = 2/sigma2*real(temp);       
                Jftildeb(k,k+K*(ell-1) ) = -2/sigma2*imag(temp);       
    end;%ell
end;%k

Jbarbf = Jfbarb';
Jtildebf = Jftildeb';

Jfb = [Jfbarb Jftildeb];
Jbf = [Jbarbf;Jtildebf];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Jbar(b)bar(b) & Jbar(b)tilde(b)  & Jtilde(b)Jtilde(b)  
% (all are KL x KL matrices)  - 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jbarbbarb = zeros(K*L,K*L);
Jbarbtildeb = zeros(K*L,K*L);
Jtildebtildeb = zeros(K*L,K*L);
for k = 1:K
    for ell = 1:L
                temp = Delayed_sig_f(:,k,ell)' *Delayed_sig_f(:,k,ell);
                Jbarbbarb(k+K*(ell-1),k+K*(ell-1)) = 2/sigma2*real(temp);       
                Jbarbtildeb(k+K*(ell-1),k+K*(ell-1)) = 0;       
    end;%ell
end;%k

Jtildebtildeb = Jbarbbarb;
Jtildebbarb =  -Jbarbtildeb;
Jbb =   [Jbarbbarb Jbarbtildeb;Jtildebbarb Jtildebtildeb];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (Known_Unknown_Sig == 0) % 1- Known Signal, 0 - Unknown Signals  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define the ds derivative 
for k = 1:K
    for ell = 1:L
        Delayed_sig_f_dS(:,k,ell) = exp(-j*2*pi*TOA(ell,k)*freq_vec');
    end;%ell
end;%k

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Calaculating Jpbar(s) (2Q x K*Ns matrix)  & Jptilde(s) (2Q x K*Ns matrix)
% - 10
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jxbars = zeros(Q,K*Ns);
Jxtildes = zeros(Q,K*Ns);
Jybars  = zeros(Q,K*Ns);
Jytildes = zeros(Q,K*Ns);

for k = 1:K
    for ell = 1:L
        temp1 = abs(b_vec(:,k,ell))^2 * temp_Diff_sig_dX(:,k,ell)' * diag(Delayed_sig_f_dS(:,k,ell));
        temp2 = abs(b_vec(:,k,ell))^2 * temp_Diff_sig_dY(:,k,ell)' * diag(Delayed_sig_f_dS(:,k,ell));
        Jxbars(1,1+(k-1)*Ns:k*Ns)  = Jxbars(1,1+(k-1)*Ns:k*Ns)+   2/sigma2*real(temp1);       
        Jxtildes(1,1+(k-1)*Ns:k*Ns) =   Jxtildes(1,1+(k-1)*Ns:k*Ns) -   2/sigma2*imag(temp1);       
        Jybars(1,1+(k-1)*Ns:k*Ns)  = Jybars(1,1+(k-1)*Ns:k*Ns) +   2/sigma2*real(temp2);       
        Jytildes(1,1+(k-1)*Ns:k*Ns) =   Jytildes(1,1+(k-1)*Ns:k*Ns) -   2/sigma2*imag(temp2);       
    end;%ell
end;%k

Jps =   [Jxbars,Jxtildes;Jybars,Jytildes];
Jsp =   Jps';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculating Jfbar(s) & Jftilde(s)  )   (K x K*Ns matrix)  - 5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jfbars = zeros(K,K*Ns);
Jftildes = zeros(K,K*Ns);

for k = 1:K
    for ell = 1:L
                temp  = abs(b_vec(:,k,ell))^2 * temp_Diff_sig_dF(:,k,ell)'*diag(Delayed_sig_f_dS(:,k,ell));
                Jfbars(k,1+(k-1)*Ns:k*Ns)  =  Jfbars(k,1+(k-1)*Ns:k*Ns) + 2/sigma2*real(temp);       
                Jftildes(k,1+(k-1)*Ns:k*Ns) =   Jftildes(k,1+(k-1)*Ns:k*Ns) -2/sigma2*imag(temp);       
    end;%ell
end;%k

Jfs =   [Jfbars Jftildes];
Jsf=Jfs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculating Jbar(b)bar(s) & Jbar(b)tilde(s)  & Jtilde(b)Jtilde(s)  (KL x K*Ns matrix)  - 7
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jbarbbars = zeros(K*L,K*Ns);
Jbarbtildes = zeros(K*L,K*Ns);
Jtildebtildes = zeros(K*L,K*Ns);
for k = 1:K
    for ell = 1:L
                temp = b_vec(:,k,ell) * Delayed_sig_f(:,k,ell)' * diag(Delayed_sig_f_dS(:,k,ell));
                Jbarbbars(k+K*(ell-1),1+(k-1)*Ns:k*Ns) =  2/sigma2*real(temp);       
                Jbarbtildes(k+K*(ell-1),1+(k-1)*Ns:k*Ns) =  -2/sigma2*imag(temp);       
    end;%ell
end;%k
Jtildebtildes = Jbarbbars;
Jtildebbars =  -Jbarbtildes;
Jbs =   [Jbarbbars Jbarbtildes;Jtildebbars Jtildebtildes];
Jsb=Jbs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculating Jbar(s)bar(s) & Jbar(s)tilde(s)  & Jtilde(s)Jtilde(s)  (K*Ns x K*Ns matrix)  - 4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jbarsbars = zeros(K*Ns,K*Ns);
Jbarstildes = zeros(K*Ns,K*Ns);
Jtildestildes = zeros(K*Ns,K*Ns);
for k = 1:K
    for ell = 1:L
                temp  = abs(b_vec(:,k,ell))^2*eye(Ns);
                Jbarsbars(1+(k-1)*Ns:k*Ns,1+(k-1)*Ns:k*Ns) = Jbarsbars(1+(k-1)*Ns:k*Ns,1+(k-1)*Ns:k*Ns) + 2/sigma2*real(temp);       
                Jbarstildes(1+(k-1)*Ns:k*Ns,1+(k-1)*Ns:k*Ns) = 0;       
    end;%ell
end;%k
Jtildestildes = Jbarsbars;
Jtildesbars  =  -Jbarstildes;
Jss =   [Jbarsbars Jbarstildes;Jtildesbars Jtildestildes];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vec1 = [];
for k = 1:K,    
    vec1 = [vec1,2+2*Q+1+(k-1)*Ns:k*Ns]; 
end;
vec = [vec1,K*Ns+vec1];
Jps = Jps(:,vec);     Jsp  =   Jps';
Jfs = Jfs(:,vec);       Jsf  =   Jfs';
Jbs =Jbs(:,vec);     Jsb  =   Jbs';
Jss = Jss(vec,vec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uknown transmitted signals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J   =   [Jpp                 Jpf         Jpb         Jps    ;
            Jfp                  Jff         Jfb         Jfs    ;
             Jbp                  Jbf         Jbb         Jbs   ;
             Jsp                  Jsf         Jsb         Jss   ;
            ];     
        
Temp = inv(J);
CRBp = Temp(1:2*Q,1:2*Q); 

else %if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Known transmitted signals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J   =   [Jpp                 Jpf         Jpb             ;
            Jfp                  Jff         Jfb             ;
             Jbp                  Jbf         Jbb            ;
            ];
        
Temp = inv(J);
CRBp = Temp(1:2*Q,1:2*Q); 
end; %if
