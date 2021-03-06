function CRBMat = CRB2(Fc,c, delayedSignal, delayedSignalDot, N, Fs, b, rTransmitter,v, rReceiverMat,SNRdB)

%Parameter to calculate: sl(time shifted s), sdotl(time shifted sdot),Al,sigmal

%Known Signals CRB

L = size(rReceiverMat,1);
rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
rDiffMat = rTransmitterMat-rReceiverMat;
d = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);

vNorm2 = v*v';
vNorm = sqrt(vNorm2);

cosTheta = rDiffMat(:,1)./d;
sinTheta = rDiffMat(:,2)./d;
cosPhi = rDiffMat*v'./(d*sqrt(vNorm2));
sinPhiSign = zeros(L,1);
for l=1:L
    tempA = [rDiffMat(l,:) 0];
    tempB = [v 0];
    tempC = cross(tempA,tempB);
    sinPhiSign(l) = sign(tempC(3));
end

    sinPhi = sqrt(1-cosPhi.^2).*sinPhiSign;


%Create sl, sDotl

sl = delayedSignal(:,1:N).';
sDotl = delayedSignalDot(:,1:N).';

sigmal = 10^(-SNRdB/20)*ones(1,L);

A = zeros(N,N,L);

for l=1:L
    miu = -1/c*v*rDiffMat(l,:)'/d(l);
    %frequency calculation
    hzDopplerShift = Fc*miu; 
    A(:,:,l) = diag(exp(2*pi*i*hzDopplerShift*(1:N)/Fs));
end
%We ignore for the moment the transmitter freq instability
C=eye(N);

Ts = 1/Fs;
Ntilde = diag(0:N-1);

dMldX = zeros(N,L);
dMldY = zeros(N,L);
dMldVx = zeros(N,L);
dMldVy = zeros(N,L);

for l=1:L
    dMldX(:,l) = (2*pi*i*Ts*Fc/c*vNorm/d(l)*sinPhi(l)*sinTheta(l))*b(l)*Ntilde*A(:,:,l)*C*sl(:,l)-1/c*cosTheta(l)*b(l)*A(:,:,l)*C*sDotl(:,l);
end

for l=1:L
    dMldY(:,l) = -(2*pi*i*Ts*Fc/c*vNorm/d(l)*sinPhi(l)*cosTheta(l))*b(l)*Ntilde*A(:,:,l)*C*sl(:,l)-1/c*sinTheta(l)*b(l)*A(:,:,l)*C*sDotl(:,l);
end

for l=1:L
    dMldVx(:,l) = -(2*pi*i*Fc/c*cosTheta(l))*Ts*b(l)*Ntilde*A(:,:,l)*C*sl(:,l);
end

for l=1:L
    dMldVy(:,l) = -(2*pi*i*Fc/c*sinTheta(l))*Ts*b(l)*Ntilde*A(:,:,l)*C*sl(:,l);
end


FIM = zeros(4); %(1-x,2-y,3-vx,4-vy];
%Calculation of the Fischer Information Matrix

for m=1:4
    switch m
        case {1}
            Mlm=dMldX;
        case{2}
            Mlm = dMldY;
        case {3}
            Mlm = dMldVx;
        case {4}
            Mlm = dMldVy;
    end %switch m

    for n=1:4
        switch n
            case {1}
                Mln=dMldX;
            case{2}
                Mln = dMldY;
            case {3}
                Mln = dMldVx;
            case {4}
                Mln = dMldVy;
        end %switch n

        for l=1:L
            FIM(m,n)=FIM(m,n)+(2/sigmal(l)^2)*real(Mlm(:,l)'*Mln(:,l));
        end %l
    end %n
end%m

CRBMat = inv(FIM);



