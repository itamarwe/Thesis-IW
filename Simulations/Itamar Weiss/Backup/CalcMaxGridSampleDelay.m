function maxGridSampleDelay = CalcMaxGridSampleDelay (gridX, gridY, rReceiverMat,C, Fs)
%function maxGridSampleDelay = CalcMaxGridSampleDelay (gridX, gridY,rReceiverMat,C, Fs)

L = size(rReceiverMat,1);
TRIALS = length(gridX);
%Maximum Delay Calculation
maxGridSampleDelay = 0;
%distSqrX = 0;
%distSqrY = 0;

for trial=1:TRIALS
for l=1:L
 %   distSqrX = (gridX(trial)-rReceiverMat(l,1))^2;
  %  distSqrY = (gridY(trial)-rReceiverMat(l,2))^2;
   % dist = sqrt(distSqrX+ distSqrY);
    %microsecDelay=dist/(C/1e6);

    rGridTransmitter = [gridX(trial) gridY(trial)];


    rGridTransmitterMat = ones(L,2)*[rGridTransmitter(1) 0;0 rGridTransmitter(2)];
    rDiffMat = rReceiverMat- rGridTransmitterMat;
    rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
    microsecDelay = rDistances/(C/1e6);
    sampleDelay = floor(microsecDelay*Fs/1e6);
    sampleMaxDelay = max(sampleDelay);
    
    if (sampleMaxDelay>maxGridSampleDelay)
        maxGridSampleDelay=sampleMaxDelay;
    end
end
end
