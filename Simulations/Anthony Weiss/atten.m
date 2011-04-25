function b_vec =   atten(att_mean,att_std,K,L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%-----------
% att_mean: mean of attenuations
% att_std   : standard deviation fo attenuations
% K : number of interception intervals
% L : number of receivers

% Outputs:
%--------------
% b_vec: attenuations (KxL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ell = 1:L
    for k = 1:K
        b_vec(:,k,ell) = (att_mean + att_std*randn(1))*exp(j*2*rand(1)*pi);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



