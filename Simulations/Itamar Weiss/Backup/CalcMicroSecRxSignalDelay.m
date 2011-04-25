function [maxDelay delay] = CalcMicroSecRxSignalDelay(rTransmitter, rReceiverMat,C)
%function [sampleMaxDelay sampleDelay] =CalcMicroSecRxSignalDelay(rTransmitter, rReceiverMat,C)
%Calculation of the Delay of the received signal from its true
%position in [uSec]
L = size(rReceiverMat,1);
rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
rDiffMat = rReceiverMat- rTransmitterMat;
rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
delay = rDistances/(C/1e6);
maxDelay = max(delay);
