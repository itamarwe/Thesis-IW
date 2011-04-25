function [TD,FD] = TD_FD_Known_Signals(r_vec_t,max_dop,time_vec,sig)

fig1 = 0;
K = size(r_vec_t,2);
L = size(r_vec_t,3);
Ns = length(time_vec);
Ts = time_vec(2)-time_vec(1);
N=2^5;
fd1 = (-2*max_dop:4*max_dop/(N-1):2*max_dop);
T0 = time_vec(end);

TD0 = zeros(K,L);
FD0 = zeros(K,L);

for k=1:K
    for ell=1:L
        %td1 = (-min([max_td(k),T0/2]):2*max_td/(N-1):min([max_td(k),T0/2]));
        td1 = 0:T0*.8/(N-1):T0*.4;
        sig1 = sig(k,:).';
        sig2 = r_vec_t(:,k,ell);
        sigf2 = repmat(sig2,[1,length(fd1)]).*exp(-1i*2*pi*time_vec'*fd1); % freq shift (in time domain)
        temp0 = fftshift(fft(sigf2),1); % freq domain
        temp1 = repmat(temp0,[1,1,length(td1)]);
        clear temp2
        temp2(:,1,:) = exp(1i*2*pi*(-Ns/2:Ns/2-1)'./Ns*(td1/Ts)); % delay in freq domain
        temp3 = repmat(temp2,[1,length(fd1),1]);
        temp4 = temp3.*temp1;

        temp5 = fftshift(fft(sig1),1);
        temp6 = repmat(temp5,[1, length(fd1), length(td1)]);
        temp7 = abs(sum(temp6.*conj(temp4),1));
        CAF = squeeze(temp7);
        [f,t] = find(CAF==max(max(CAF)));

        TD0(k,ell) = td1(t(1)) ;
        FD0(k,ell) = fd1(f(1)) ;

        fd2 = FD0(k,ell) + (-(fd1(2)-fd1(1))/2 : (fd1(2)-fd1(1))/(2*N): (fd1(2)-fd1(1))/2 ) ;
        td2 = TD0(k,ell) + (-(td1(2)-td1(1))/2 : (td1(2)-td1(1))/(2*N): (td1(2)-td1(1))/2 ) ;

        sigf2 = repmat(sig2,[1,length(fd2)]).*exp(-1i*2*pi*time_vec'*fd2); % freq shift (in time domain)
        temp0 = fftshift(fft(sigf2),1); % freq domain
        temp1 = repmat(temp0,[1,1,length(td2)]);
        clear temp2 temp3 temp4 temp6 temp7 CAF
        temp2(:,1,:) = exp(1i*2*pi*(-Ns/2:Ns/2-1)'./Ns*(td2/Ts)); % delay in freq domain
        temp3 = repmat(temp2,[1,length(fd2),1]);
        temp4 = temp3.*temp1;

        %temp5 = fftshift(fft(sig1),1);
        temp6 = repmat(temp5,[1, length(fd2), length(td2)]);
        temp7 = abs(sum(temp6.*conj(temp4),1));
        CAF = squeeze(temp7);
        [f,t] = find(CAF==max(max(CAF)));

        TD0(k,ell) = td2(t(1)) ;
        FD0(k,ell) = fd2(f(1)) ;

        fd3 = FD0(k,ell) + (-(fd2(2)-fd2(1))/2 : (fd2(2)-fd2(1))/(2*N): (fd2(2)-fd2(1))/2 ) ;
        td3 = TD0(k,ell) + (-(td2(2)-td2(1))/2 : (td2(2)-td2(1))/(2*N): (td2(2)-td2(1))/2 ) ;

        sigf2 = repmat(sig2,[1,length(fd3)]).*exp(-1i*2*pi*time_vec'*fd3); % freq shift (in time domain)
        temp0 = fftshift(fft(sigf2),1); % freq domain
        temp1 = repmat(temp0,[1,1,length(td3)]);
        clear temp2 temp3 temp4 temp6 temp7 CAF
        temp2(:,1,:) = exp(1i*2*pi*(-Ns/2:Ns/2-1)'./Ns*(td3/Ts)); % delay in freq domain
        temp3 = repmat(temp2,[1,length(fd3),1]);
        temp4 = temp3.*temp1;

        %temp5 = fftshift(fft(sig1),1);
        temp6 = repmat(temp5,[1, length(fd3), length(td3)]);
        temp7 = abs(sum(temp6.*conj(temp4),1));
        CAF = squeeze(temp7);
        [f,t] = find(CAF==max(max(CAF)));


        if f(1)~=1 && f(1)~=length(fd3)
            [peakf] = peak_identification(fd3, CAF(:,t(1)), f(1));
            FD0(k,ell) = peakf; %fd3(f(1)) ;
        else
            FD0(k,ell) = fd3(f(1)) ;
        end
        if t(1)~=1 && t(1)~=length(td3)
            [peakt] = peak_identification(td3, CAF(f(1),:), t(1));
            TD0(k,ell) = peakt; %td3(t(1)) ;
        else
            TD0(k,ell) = td3(t(1)) ;
        end

        if fig1==1 
            figure
            mesh(CAF)
            shading interp
            figure
            contour(CAF)
            keyboard
        end
        clear f t
    end
end
TD=TD0;
FD=FD0;
%TD0,FD0,

            
