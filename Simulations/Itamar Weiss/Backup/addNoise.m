function [receivedSignal] = addNoise(cleanReceivedSignal,SNRdB);
[L,N] = size(cleanReceivedSignal);
noise = zeros(L,N);
%Receiver Noise
for l=1:L
    a = sqrt(cleanReceivedSignal(l,:)*cleanReceivedSignal(l,:)');
    b = (10^(-SNRdB/10));
    noise(l,:) = a*b*(1/sqrt(2))*(randn(1,N)+i*randn(1,N))/sqrt(N);
    receivedSignal = cleanReceivedSignal+noise;
end