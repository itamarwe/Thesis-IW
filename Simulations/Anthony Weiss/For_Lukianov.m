dLR_DPD =  2000;
Xsearch_DPD = [-dLR_DPD,-dLR_DPD/2,0,dLR_DPD/2,dLR_DPD]+EmitterPos(1)+(rand-.5)*10;
Ysearch_DPD = [-dLR_DPD,-dLR_DPD/2,0,dLR_DPD/2,dLR_DPD]+EmitterPos(2)+(rand-.5)*10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIRECT EMITTER LOCALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter = 0;     
while dLR_DPD > .1
    counter = counter+1;
    [EmitPosEst_DPD,Xest_DPD,Yest_DPD,CF_DPD]  = DPD_Known_Signals(time_vec,Xsearch_DPD,Ysearch_DPD,...
        ReceiverPos,ReceiverVel,K,L,Q,fc,c,r_vec_t,Known_Unknown_Sig,CF1,sig);
    dLR_DPD=dLR_DPD*.8;
    clear Xsearch_DPD Ysearch_DPD
    Xsearch_DPD = [-dLR_DPD,0,dLR_DPD]+Xest_DPD;
    Ysearch_DPD = [-dLR_DPD,0,dLR_DPD]+Yest_DPD;
end