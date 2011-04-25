function [EmitPosEst,Xest,Yest] = DPD_Known_Signals(time_vec,Xsearch,Ysearch,ReceiverPos,ReceiverVel,...
    K,L,Q,fc,c,r_vec_t0,CF1,sig)

CF = zeros(length(Xsearch),length(Ysearch));
Ts=time_vec(2)-time_vec(1);
Ns=length(time_vec);
%Vk=zeros(Ns,L);
for x_ind = 1:length(Xsearch)
       for y_ind = 1:length(Ysearch)
                EmitPostmp = [Xsearch(x_ind),Ysearch(y_ind)];
                Dop_shift = Create_g_shift(EmitPostmp,ReceiverPos,ReceiverVel,fc,c,L,K,Q);
                %CF_temp = 0;
                temp1 = zeros(L,K);
                for k = 1:K
                    TOA_rel=zeros(1,L);
                    for ell = 1:L
                        TOA_rel(ell) = norm(ReceiverPos(k,:,ell) - EmitPostmp)/c; %TOA_rel = Relative distance /c
                        DD = Dop_shift(ell,k,1);
                        sig_t = r_vec_t0(:,k,ell).*exp(-1i*2*pi*DD*time_vec');
                        DT = TOA_rel(ell);
                        temp = fftshift(fft(sig_t)).*exp(1i*2*pi*(-Ns/2:Ns/2-1)'./Ns*(DT/Ts));
                        %Vk(:,ell) = temp;
                        temp1(ell,k) = abs(conj(fftshift(fft(sig(k,:))))*temp)^2;
                    end;%ell
                    %Lambda  = svd(Vk'*Vk); % (Eq.22 - old)
                    %CF_temp =  CF_temp + abs(Lambda(1));
                    
                end;%k
                %CF(x_ind,y_ind) = CF_temp;
                CF(x_ind,y_ind) = sum(sum(temp1));
       end; % y_ind
end; % x_ind
if CF1==1
 figure
 surf(CF), shading interp
  title('DPD CF')
 figure
 contour(CF), 
  title('DPD CF')
 keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the optimum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u1,v1] = find(CF == max(max(CF))); 
Xest = Xsearch(u1(1)); 
Yest = Ysearch(v1(1));
EmitPosEst = [Xest,Yest]';
% figure,surf(CF), shading interp

