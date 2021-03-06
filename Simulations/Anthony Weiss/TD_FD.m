function [TD,FD] = TD_FD(r_vec_t,max_dop,max_td,time_vec)

fig1 = 0;
K = size(r_vec_t,2);
Ns = length(time_vec);
Ts = time_vec(2)-time_vec(1);
N=2^5;
fd1 = (-2*max_dop:4*max_dop/(N-1):2*max_dop);
T0 = time_vec(end);

TD = zeros(K,1);
FD = zeros(K,1);

for k=1:K
    td1 = (-min([max_td(k),T0/2]):2*max_td/(N-1):min([max_td(k),T0/2]));
    sig1 = r_vec_t(:,k,1);
    sig2 = r_vec_t(:,k,2);
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
    
    TD(k,:) = td1(t(1)) ;
    FD(k,:) = fd1(f(1)) ;
    
    fd2 = FD(k) + (-(fd1(2)-fd1(1))/2 : (fd1(2)-fd1(1))/(2*N): (fd1(2)-fd1(1))/2 ) ;
    td2 = TD(k) + (-(td1(2)-td1(1))/2 : (td1(2)-td1(1))/(2*N): (td1(2)-td1(1))/2 ) ;
    
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
    
    TD(k,:) = td2(t(1)) ;
    FD(k,:) = fd2(f(1)) ;
    
    fd3 = FD(k) + (-(fd2(2)-fd2(1))/2 : (fd2(2)-fd2(1))/(2*N): (fd2(2)-fd2(1))/2 ) ;
    td3 = TD(k) + (-(td2(2)-td2(1))/2 : (td2(2)-td2(1))/(2*N): (td2(2)-td2(1))/2 ) ;
    
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
        FD(k,:) = peakf; %fd3(f(1)) ;
    else
        FD(k,:) = fd3(f(1)) ;
    end
    if t(1)~=1 && t(1)~=length(td3)
        [peakt] = peak_identification(td3, CAF(f(1),:), t(1));
        TD(k,:) = peakt; %td3(t(1)) ;
    else
        TD(k,:) = td3(t(1)) ;
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
%TD,FD,

            
