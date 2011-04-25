function [EmitPosEst_ST,Xest_ST,Yest_ST]  =  STPosEst0(FD,TD,XsearchST,YsearchST,K,L,ReceiverPos,ReceiverVel,fc,c,CF1,CRB_tau,CRB_f)       


w_toa = 1/CRB_tau;  % the weight of the toa data (1/sigma_toa)
w_dop = 1/CRB_f; % the weight of the dop data (1/sigma_dop)

CF_ST = zeros(length(XsearchST),length(YsearchST));


for x_ind  =   1:length(XsearchST)
      for y_ind  =   1:length(YsearchST)
          
          TOA_tmp = zeros(L,K);
          EmitterRelPos_tmp = zeros(2,L,K);
          Vrel_tmp = zeros(L,K);
          Doppler_shift_tmp = zeros(L,K);
          EmitPostmp  =   [XsearchST(x_ind),YsearchST(y_ind)];
          
          for k = 1:K
              for ell = 1:L
                  TOA_tmp(ell,k) = norm(   ReceiverPos(k,:,ell)  -   EmitPostmp  ) /c;
                  EmitterRelPos_tmp(:,ell,k) = [EmitPostmp(1) - ReceiverPos(k,1,ell), EmitPostmp(2) - ReceiverPos(k,2,ell)];
                  Vrel_tmp(ell,k) = ReceiverVel(k,:,ell)*(EmitterRelPos_tmp(:,ell,k)./norm(EmitterRelPos_tmp(:,ell,k),2));
                  Doppler_shift_tmp(ell,k) = fc*Vrel_tmp(ell,k)/c; %  (Eq.3- new)
              end%ell;
          end%k;
          CF_ST(x_ind,y_ind) = w_toa*sum(sum((TOA_tmp - TD').^2)) + w_dop*sum(sum(( Doppler_shift_tmp - FD'  ).^2))  ;

      end; %x_ind
end; %y_ind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the optimum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u,v]     =   find(CF_ST == min(min(CF_ST))); 
Xest_ST   =   mean(XsearchST(u)); 
Yest_ST   =   mean(YsearchST(v));

EmitPosEst_ST  =   [Xest_ST,Yest_ST]';
if CF1==1
    figure
    surf(CF_ST), shading interp
    title('STEIN CF')
    figure
    contour(CF_ST), 
     title('STEIN CF')
    keyboard
end
