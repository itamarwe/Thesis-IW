function varargout = DOPGUI(varargin)
% DOPGUI M-file for DOPGUI.fig
%      DOPGUI, by itself, creates a new DOPGUI or raises the existing
%      singleton*.
%
%      H = DOPGUI returns the handle to a new DOPGUI or the handle to
%      the existing singleton*.
%
%      DOPGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DOPGUI.M with the given input arguments.
%
%      DOPGUI('Property','Value',...) creates a new DOPGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DOPGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DOPGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DOPGUI

% Last Modified by GUIDE v2.5 24-Apr-2010 19:55:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DOPGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DOPGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DOPGUI is made visible.
function DOPGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DOPGUI (see VARARGIN)

% Choose default command line output for DOPGUI
handles.output = hObject;

handles.N = str2num(get(handles.editN,'String')); %Number of samples used for estimation
handles.L = str2num(get(handles.editL,'String')); %Number of receivers
handles.Fs = str2num(get(handles.editFs,'String')); %Sampling Frequency [Hz]
handles.Trials = str2num(get(handles.editTrials,'String')); %Number of Monte Carlo Trials
handles.R  = str2num(get(handles.editR,'String')) %In Circular Receiver Array - The radius of the circle
handles.Px = str2num(get(handles.editPx,'String')); %The transmitter's X position
handles.Py = str2num(get(handles.editPy,'String')); %The transmitter's Y position
handles.Vx = str2num(get(handles.editVx,'String')); %The transmitter's X velocity
handles.Vy = str2num(get(handles.editFVy,'String')); %The transmitter's Y velocity
handles.Fc = str2num(get(handles.editFc,'String')); %The transmitter's Carrier Frequency
handles.Fsig = str2num(get(handles.editFsig,'String')); %The signal's modulation frequency

%Additional Parameters:
%Signal's Bandwidth
%Circular/Linear Array
%Grid Parameters
%Known/Unknown signals
%Gradient Search XY
%Gradient Search Vx Vy
%SNR

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DOPGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DOPGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editFVx_Callback(hObject, eventdata, handles)
% hObject    handle to editFVx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFVx as text
%        str2double(get(hObject,'String')) returns contents of editFVx as a double


% --- Executes during object creation, after setting all properties.
function editFVx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFVx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to editFVy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFVy as text
%        str2double(get(hObject,'String')) returns contents of editFVy as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFVy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFPx_Callback(hObject, eventdata, handles)
% hObject    handle to editFPx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFPx as text
%        str2double(get(hObject,'String')) returns contents of editFPx as a double


% --- Executes during object creation, after setting all properties.
function editFPx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFPx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFPy_Callback(hObject, eventdata, handles)
% hObject    handle to editFPy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFPy as text
%        str2double(get(hObject,'String')) returns contents of editFPy as a double


% --- Executes during object creation, after setting all properties.
function editFPy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFPy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPx_Callback(hObject, eventdata, handles)
% hObject    handle to editPx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPx as text
%        str2double(get(hObject,'String')) returns contents of editPx as a double


% --- Executes during object creation, after setting all properties.
function editPx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPy_Callback(hObject, eventdata, handles)
% hObject    handle to editPy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPy as text
%        str2double(get(hObject,'String')) returns contents of editPy as a double


% --- Executes during object creation, after setting all properties.
function editPy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editVx_Callback(hObject, eventdata, handles)
% hObject    handle to editVx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVx as text
%        str2double(get(hObject,'String')) returns contents of editVx as a double


% --- Executes during object creation, after setting all properties.
function editVx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFVy_Callback(hObject, eventdata, handles)
% hObject    handle to editFVy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFVy as text
%        str2double(get(hObject,'String')) returns contents of editFVy as a double


% --- Executes during object creation, after setting all properties.
function editFVy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFVy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTrials_Callback(hObject, eventdata, handles)
% hObject    handle to editTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTrials as text
%        str2double(get(hObject,'String')) returns contents of editTrials as a double


% --- Executes during object creation, after setting all properties.
function editTrials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFs_Callback(hObject, eventdata, handles)
% hObject    handle to editFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFs as text
%        str2double(get(hObject,'String')) returns contents of editFs as a double


% --- Executes during object creation, after setting all properties.
function editFs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editN_Callback(hObject, eventdata, handles)
% hObject    handle to editN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editN as text
%        str2double(get(hObject,'String')) returns contents of editN as a double


% --- Executes during object creation, after setting all properties.
function editN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFc_Callback(hObject, eventdata, handles)
% hObject    handle to editFc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFc as text
%        str2double(get(hObject,'String')) returns contents of editFc as a double


% --- Executes during object creation, after setting all properties.
function editFc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editL_Callback(hObject, eventdata, handles)
% hObject    handle to editL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editL as text
%        str2double(get(hObject,'String')) returns contents of editL as a double


% --- Executes during object creation, after setting all properties.
function editL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editR_Callback(hObject, eventdata, handles)
% hObject    handle to editR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editR as text
%        str2double(get(hObject,'String')) returns contents of editR as a double


% --- Executes during object creation, after setting all properties.
function editR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function run(handles)
%Simulation 02 - Model

%*****************
%*** init ***
%*****************
C = 3e08; %Speed of Signal Propagation (Usually speed of light)[m/s]
N = str2num(get(handles.editN,'String')); %Number of samples used for estimation
L = str2num(get(handles.editL,'String')); %Number of receivers
Fs = str2num(get(handles.editFs,'String')); %Sampling Frequency [Hz]
TRIALS = str2num(get(handles.editTrials,'String')); %Number of Monte Carlo Trials


%*** Receivers Positions ***

%Spread circulary around the axis center
R = str2num(get(handles.editR,'String'));% [m] Distance From the axis center

rReceiverMat = R*[sin(((0:(L-1))*2*pi/L))' cos(((0:(L-1))*2*pi/L))' ]; %Positions of the receivers in the form [x1 y1; x2 y2;...etc]

%Transmitter Position, Velocity, Nominal Frequency, Bandwidth
rTransmitter = [str2num(get(handles.editPx,'String')) str2num(get(handles.editPy,'String'))]; %Transmitter position [x,y] [m]
v = [str2num(get(handles.editVx,'String')) str2num(get(handles.editVy,'String'))];%Transmitter velocity [vx, vy] [m/s]
Fc = str2num(get(handles.editFc,'String')); % Carrier Frequency[Hz]

%for l=1:L
%    for n=1:N
%       vA = exp(2*pi*i*hzDopplerShift(l)*(1:N)/Fs)';
%       receivedSignal(l,:) = b(l)*(vA*delayedSignal(l,:).').'+noise(l,:);
%end

%    end
%end


%/////////////////////////////////////
%*************************************
%** Estimation ***
%*************************************
%/////////////////////////////////////


%************************************
% GRID SEARCH
%************************************

%Grid Search Init


%Determine Grid According to the user's preference
if (get(handles.radiobuttonGridGeneral,'Value'));
gridX = rTransmitter(1) + [randn(1,TRIALS-1)*100 0];
gridY = rTransmitter(2) + [randn(1,TRIALS-1)*100 0];
% gridX = zeros(1,TRIALS)+100;
% gridY = zeros(1,TRIALS)+150;
gridVx = v(1)*([1+randn(1,TRIALS-1)*1 1]);
gridVy = v(2)*([1+randn(1,TRIALS-1)*1 1]);
% gridVx = v(1)*ones(1,TRIALS)+10;
% gridVy = v(2)*ones(1,TRIALS)+12;
elseif (get(handles.radiobuttonGridFixedPosition,'Value'));
gridX = str2num(get(handles.editFPx,'String'))*ones(1,TRIALS);
gridY = str2num(get(handles.editFPy,'String'))*ones(1,TRIALS);
% gridX = zeros(1,TRIALS)+100;
% gridY = zeros(1,TRIALS)+150;
gridVx = v(1)*([1+randn(1,TRIALS-1)*0.1 1]);
gridVy = v(2)*([1+randn(1,TRIALS-1)*0.1 1]);
% gridVx = v(1)*ones(1,TRIALS)+10;
% gridVy = v(2)*ones(1,TRIALS)+12;
elseif (get(handles.radiobuttonGridFixedVelocity,'Value'));
gridX = rTransmitter(1) + [randn(1,TRIALS-1)*100 0];
gridY = rTransmitter(2) + [randn(1,TRIALS-1)*100 0];

gridVx = str2num(get(handles.editFVx,'String'))*ones(1,TRIALS);
gridVy = str2num(get(handles.editFVy,'String'))*ones(1,TRIALS);
end


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

%   maxGridSampleDelay = floor(maxGridMicrosecDelay*Fs/1e6);

%//////////////////////////////////
%***********************************
%*** Received Downshifted Signal ***
%***********************************
%//////////////////////////////////

%Calculation of the Longest Delay
 rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
 rDiffMat = rReceiverMat- rTransmitterMat;
 rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
 microsecDelay = rDistances/(C/1e6);
 sampleDelay = floor(microsecDelay*Fs/1e6);
 sampleMaxDelay = max(sampleDelay);



%Original Signal
Fsig = str2num(get(handles.editFsig,'String'));
%s = randn(1,N+sampleMaxDelay).*exp(j*randn(1,N+sampleMaxDelay));
s = 100*sin(2*pi*Fsig*(1:(N+maxGridSampleDelay+sampleMaxDelay))/Fs);

%Complex Path Attenuation
b = randn(1,L)+ randn(1,L)*i;
%b=ones(1,L);
%Receiver Noise
noise = 0*randn(L,N+maxGridSampleDelay);

%Delayed Signal
delayedSignal = zeros(L,N+maxGridSampleDelay);

for j=1:L
    delayedSignal(j,:) = s((1+sampleMaxDelay-sampleDelay(j)):(N+maxGridSampleDelay+sampleMaxDelay-sampleDelay(j)));
end

%A = zeros(N,N,L);

miu = zeros(L,1);
hzDopplerShift = zeros(L,1);
%Doppler Shift
for l=1:L
    %miu calculation
    miu(l) = -1/C*v*rDiffMat(l,:)'/sqrt(rDiffMat(l,:)*rDiffMat(l,:)');
    %frequency calculation
    hzDopplerShift(l) = Fc*miu(l);
    %Al Matrix
    %A(:,:,l) = diag(exp(2*pi*i*hzDopplerShift(l)*(1:N)/Fs));
end


%Received Signal
receivedSignal = zeros(L,N+maxGridSampleDelay);

for l=1:L
    vA=exp(2*pi*i*hzDopplerShift(l)*(1:(N+maxGridSampleDelay))/Fs);
    receivedSignal(l,:) = b(l)*vA.*delayedSignal(l,:)+noise(l,:);
end



lambdaMax = zeros(1,TRIALS);
waitbarHandle = waitbar(0,'Processing');
for trial=1:TRIALS

    if mod(trial,10)==0
        waitbar(trial/TRIALS,waitbarHandle);
    end

    rTransmitter = [gridX(trial) gridY(trial)];
    v = [gridVx(trial) gridVy(trial)];


    rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
    rDiffMat = rReceiverMat- rTransmitterMat;
    rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
    microsecDelay = rDistances/(C/1e6);
    sampleDelay = floor(microsecDelay*Fs/1e6);
    sampleMaxDelay = max(sampleDelay);
    %Evaluate hermitian A
    %    hermitianA = zeros(N,N,L);
    miu = zeros(L,1);
    dopplerFixedSignal = zeros(L,N+maxGridSampleDelay);
    for l=1:L
        %miu calculation
        miu(l) = -1/C*v*rDiffMat(l,:)'/rDistances(l);
        %frequency calculation
        hzDopplerShift = Fc*miu(l);
        %Al Matrix
        %hermitianA(:,:,l) = diag(exp(2*pi*i*hzDopplerShift*(1:N)/Fs))';
        %dopplerFixedSignal(l,:) = (hermitianA(:,:,l)*receivedSignal(l,:).').';
        vHermitianA = exp(-2*pi*i*hzDopplerShift*(1:N+maxGridSampleDelay)/Fs);
        dopplerFixedSignal(l,:) =(vHermitianA.*receivedSignal(l,:)).';
    end

   % Ntilda = N-sampleMaxDelay;
    %Evaluate V
    V = zeros(L,N);
    for l=1:L
        V(l,:) = dopplerFixedSignal(l,(1+sampleDelay(l)):(N+sampleDelay(l)));
    end

    %Evaluate Qtilda


    Qtilda = V*V';
    %Evaluate LambdaMax(Qtilda)
    lambdaMax(trial) = eigs(Qtilda,1);
end

%*** Plots ****

% %Transmitter And Receivers Position
% figure;
% hold on;
% plot(rTransmitter(1),rTransmitter(2),'+r');
% plot(rReceiverMat(:,1),rReceiverMat(:,2),'*b');
% hold off;
% grid on;
% legend('Transmitter Position','Receivers Position');
% 
% %Frequency Domain plot of the signals
% figure;
% for l=1:L
%     [Pxx(l,:),F(l,:)] = pwelch(receivedSignal(l,:),[],[],[],Fs);
% end
% plot(F',10*log10(Pxx)');
% title('Frequency Domain plot of the signals');

%figure;
if ~(get(handles.radiobuttonGridFixedPosition,'Value'))
[nGridX nGridY] = meshgrid(linspace(min(gridX),max(gridX),200),linspace(min(gridY),max(gridY),200));
nLambdaMax = griddata(gridX, gridY, lambdaMax,nGridX, nGridY);
surfc(handles.axesP,nGridX,nGridY,nLambdaMax)
title(handles.axesP,'Position Cost Function');
else plot(handles.axesP,1,1);
end

%figure
if ~(get(handles.radiobuttonGridFixedVelocity,'Value'))
[nGridVx nGridVy] = meshgrid(linspace(min(gridVx),max(gridVx),200),linspace(min(gridVy),max(gridVy),200));
nLambdaMax = griddata(gridVx, gridVy, lambdaMax,nGridVx, nGridVy);
surfc(handles.axesV,nGridVx,nGridVy,nLambdaMax);
title(handles.axesV,'Velocity Cost Function');
else plot(handles.axesV,1,1);
end

drawnow;
close(waitbarHandle);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run(handles);

% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function editFsig_Callback(hObject, eventdata, handles)
% hObject    handle to editFsig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFsig as text
%        str2double(get(hObject,'String')) returns contents of editFsig as a double


% --- Executes during object creation, after setting all properties.
function editFsig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFsig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editVy_Callback(hObject, eventdata, handles)
% hObject    handle to editVy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVy as text
%        str2double(get(hObject,'String')) returns contents of editVy as a double


% --- Executes during object creation, after setting all properties.
function editVy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


