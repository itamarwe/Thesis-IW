function maxGridSampleDelay = CalcMaxGridSampleDelay2 (cX,cY, gridResolution, gridPoints, rReceiverMat,C, Fs)
%function maxGridSampleDelay = CalcMaxGridSampleDelay2 (cX,cY,gridResolution, gridPoints, rReceiverMat,C, Fs)

L = size(rReceiverMat,1);

gridX = ((-gridPoints/2):1:(gridPoints/2))*gridResolution(1)+cX;
gridY = ((-gridPoints/2):1:(gridPoints/2))*gridResolution(2)+cY;

%Maximum Delay Calculation
maxGridSampleDelay = 0;


for gX = gridX
    for gY = gridY
        for l=1:L
            %   distSqrX = (gridX(trial)-rReceiverMat(l,1))^2;
            %  distSqrY = (gridY(trial)-rReceiverMat(l,2))^2;
            % dist = sqrt(distSqrX+ distSqrY);
            %microsecDelay=dist/(C/1e6);

            rGridTransmitter = [gX gY];


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
end
