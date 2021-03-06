function [gridX, gridY, gridVx, gridVy] = DefineGrid(rTransmitter, v, Trials, method, FPxOrVx, FPyOrVy)
%Method:1 - GeneralSearch 2 - Fixed Position 3 - Fixed Velocity
%Determine Grid According to the user's preference
if (method==1)
    %Grid search in the position and velocity subspaces
    gridX = rTransmitter(1) + [randn(1,Trials-1)*100 0];
    gridY = rTransmitter(2) + [randn(1,Trials-1)*100 0];
    % gridX = zeros(1,TRIALS)+100;
    % gridY = zeros(1,TRIALS)+150;
    gridVx = v(1)+([1+randn(1,Trials-1)*200 0]);
    gridVy = v(2)+([1+randn(1,Trials-1)*200 0]);
    % gridVx = v(1)*ones(1,TRIALS)+10;
    % gridVy = v(2)*ones(1,TRIALS)+12;
elseif (method==2)  
    %Grid search only in the velocity subspace, with fixed position
    gridX = FPxOrVx*ones(1,Trials);
    gridY = FPyOrVy*ones(1,Trials);
    % gridX = zeros(1,TRIALS)+100;
    % gridY = zeros(1,TRIALS)+150;
    gridVx = v(1)+[(rand(1,Trials-1)-0.5)*200 0];
    gridVy = v(2)+[(rand(1,Trials-1)-0.5)*200 0];
    % gridVx = v(1)*ones(1,TRIALS)+10;
    % gridVy = v(2)*ones(1,TRIALS)+12;
elseif (method==3)
    %Grid search only in the position subspace, with fixed velocity
    gridX = rTransmitter(1) + [(rand(1,Trials-1)-0.5)*500 0];
    gridY = rTransmitter(2) + [(rand(1,Trials-1)-0.5)*500 0];

    gridVx = FPxOrVx*ones(1,Trials);
    gridVy = FPyOrVy*ones(1,Trials);
end

