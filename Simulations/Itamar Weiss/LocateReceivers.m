function rReceiverMat = LocateReceivers(R,L, method)
    %Locating the L receivers.
    %The center of the array is in [0,0] with distance parameter R[m]
    %Method: 1 - Circular Array with radius R
    %Method 2 - Linear Array with length 2R all elements are located on the
    %X axis
    %rReceiverMat - Positions of the receivers in the form [x1 y1; x2 y2;...etc]
    if (method==1)  
        rReceiverMat = R*[sin(((0:(L-1))*2*pi/L))' cos(((0:(L-1))*2*pi/L))' ]; %Positions of the receivers in the form [x1 y1; x2 y2;...etc]
    end
    
    if (method==2)  
        rReceiverMat = [R*linspace(-1,1,L);zeros(1,L)]';
    end
    
