function [sampleMaxDelay sampleDelay] = CalcRxSignalDelay(rTransmitter, rReceiverMat,C,Fs)
%function [sampleMaxDelay sampleDelay] = CalcMaxRxSignalDelay(rTransmitter,rReceiverMat,C,Fs)
%Calculation of the Delay of the received signal from its true
%position
L = size(rReceiverMat,1);
rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
rDiffMat = rReceiverMat- rTransmitterMat;
rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
microsecDelay = rDistances/(C/1e6);
sampleDelay = floor(microsecDelay*Fs/1e6);
sampleMaxDelay = max(sampleDelay);
