function varargout = measGUI(varargin)
% MEASGUI MATLAB code for measGUI.fig
%      MEASGUI, by itself, creates a new MEASGUI or raises the existing
%      singleton*.
%
%      H = MEASGUI returns the handle to a new MEASGUI or the handle to
%      the existing singleton*.
%
%      MEASGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEASGUI.M with the given input arguments.
%
%      MEASGUI('Property','Value',...) creates a new MEASGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before measGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to measGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help measGUI

% Last Modified by GUIDE v2.5 02-Jun-2016 14:35:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @measGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @measGUI_OutputFcn, ...
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


% --- Executes just before measGUI is made visible.
function measGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to measGUI (see VARARGIN)

% Choose default command line output for measGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% default date is current date
today = datestr(date, 29);
dashes = strfind(today, '-');
currYear = today(3:4);
currMonth = num2str(str2num(today(dashes(1)+1:dashes(2)-1)));
currDay = num2str(str2num(today(dashes(2)+1:end)));

% check if settings file exists and load it, otherwise start off all values
% as default
settingsFileName = 'C:\Users\Livia\Desktop\IVF\Code\EmbryoProject\GUIsettings.mat';
if exist(settingsFileName, 'file')
    
    % load previous state
    load(settingsFileName);
    setappdata(handles.GUI, 'currDate', currDate);
    setappdata(handles.GUI, 'dateU', dateU);
    setappdata(handles.GUI, 'measType', measType);
    setappdata(handles.GUI, 'pipRefExists', pipRefExists);
    setappdata(handles.GUI, 'embryoNum', embryoNum);
    setappdata(handles.GUI, 'pipSize', pipSize);
    setappdata(handles.GUI, 'manualMeasure', manualMeasure);
    setappdata(handles.GUI, 'filePathRaw', filePathRaw);
    setappdata(handles.GUI, 'filePathProc', filePathProc);
    setappdata(handles.GUI, 'videoLoaded', videoLoaded);
    setappdata(handles.GUI, 'frames', []);
    setappdata(handles.GUI, 'lastFrame', lastFrame);
    setappdata(handles.GUI, 'currFrame', currFrame);
    setappdata(handles.GUI, 't', t);
    setappdata(handles.GUI, 'params', params);
    setappdata(handles.GUI, 'paramsFit', paramsFit);
    setappdata(handles.GUI, 'displayFigNum', displayFigNum);
    setappdata(handles.GUI, 'plotFigNum', plotFigNum);
    setappdata(handles.GUI, 'extraFigNum', extraFigNum);
    setappdata(handles.GUI, 'zpEnter', zpEnter);
    setappdata(handles.GUI, 'cellEnter', cellEnter);
    setappdata(handles.GUI, 'alreadyCropped', alreadyCropped);
    
    % init GUI to previous state
    set(handles.DateTextEdit, 'String', getappdata(handles.GUI, 'currDate'));
    set(handles.RawDataLabel, 'String', filePathRaw);
    set(handles.ProcDataLabel, 'String', filePathProc);
    set(handles.PipSizeEdit, 'String', num2str(pipSize));
    set(handles.EmbryoNameEdit, 'String', embryoNum);

    
    
else
    setappdata(handles.GUI, 'currDate', [currMonth '-' currDay '-' currYear]);
    setappdata(handles.GUI, 'dateU', [currMonth '_' currDay '_' currYear]);
    setappdata(handles.GUI, 'measType', 'Mouse Embryo');
    setappdata(handles.GUI, 'pipRefExists', 0);
    setappdata(handles.GUI, 'embryoNum', '1');
    setappdata(handles.GUI, 'pipSize', 128);
    setappdata(handles.GUI, 'manualMeasure', 0);
    setappdata(handles.GUI, 'filePathRaw', '');
    setappdata(handles.GUI, 'filePathProc', '');
    setappdata(handles.GUI, 'videoLoaded', 0);
    setappdata(handles.GUI, 'frames', []);
    setappdata(handles.GUI, 'lastFrame', 0);
    setappdata(handles.GUI, 'currFrame', 0);
    setappdata(handles.GUI, 't', []);
    setappdata(handles.GUI, 'params', []);
    setappdata(handles.GUI, 'paramsFit', []);
    setappdata(handles.GUI, 'displayFigNum', 1);
    setappdata(handles.GUI, 'plotFigNum', 2);
    setappdata(handles.GUI, 'extraFigNum', 3);
    setappdata(handles.GUI, 'zpEnter', 0);
    setappdata(handles.GUI, 'cellEnter', 0);
    setappdata(handles.GUI, 'alreadyCropped', 0);
end

% UIWAIT makes measGUI wait for user response (see UIRESUME)
% uiwait(handles.GUI);


% --- Outputs from this function are returned to the command line.
function varargout = measGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes when user attempts to close GUI.
function GUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% write current data to a temp file to save for next time the function is
% opened
settingsFileName = 'C:\Users\Livia\Desktop\IVF\Code\EmbryoProject\GUIsettings.mat';

if isstruct(handles)
    settingsToSave = getappdata(handles.GUI)';
    
    % clear structs and potentially large variables
    if isfield(settingsToSave, 'GUIDEOptions')
        settingsToSave = rmfield(settingsToSave,'GUIDEOptions');
    end
    if isfield(settingsToSave, 'UsedByGUIData_m')
        settingsToSave = rmfield(settingsToSave,'UsedByGUIData_m');
    end
    if isfield(settingsToSave, 'frames')
        settingsToSave = rmfield(settingsToSave,'frames');
    end
    save(settingsFileName, '-struct', 'settingsToSave');
end

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on selection change in TypeSelector.
function TypeSelector_Callback(hObject, eventdata, handles)
% hObject    handle to TypeSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TypeSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TypeSelector

contents = cellstr(get(hObject,'String'));
measType = contents{get(hObject,'Value')};
setappdata(handles.GUI, 'measType', measType);



% --- Executes during object creation, after setting all properties.
function TypeSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TypeSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DateTextEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DateTextEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DateTextEdit as text
%        str2double(get(hObject,'String')) returns contents of DateTextEdit as a double

% update currDate and dateU
currDate = get(hObject, 'String');
dateU = currDate;
dateUI = strfind(currDate, '-');
dateU(dateUI) = '_';
setappdata(handles.GUI, 'currDate', currDate);
setappdata(handles.GUI, 'dateU', dateU);
measType = getappdata(handles.GUI, 'measType');

filePathRaw = ['C:\Users\Livia\Desktop\IVF\Raw Data\Videos\' measType ...
    '\videos ' currDate];
filePathProc = ['C:\Users\Livia\Desktop\IVF\Processed Data\' measType ...
    '\' currDate ' analysis'];

if ~exist(filePathRaw, 'dir')
    setappdata(handles.GUI, 'filePathRaw', '');
    set(handles.RawDataLabel, 'String', '');
    setappdata(handles.GUI, 'filePathProc', '');
    set(handles.ProcDataLabel, 'String', '');
    errordlg('No raw or processed data found for this location!');
elseif ~exist(filePathProc, 'dir')
    setappdata(handles.GUI, 'filePathRaw', filePathRaw);
    set(handles.RawDataLabel, 'String', filePathRaw);
    mkdir(filePathProc);
    setappdata(handles.GUI, 'filePathProc', filePathProc);
    set(handles.ProcDataLabel, 'String', filePathProc);
else
    setappdata(handles.GUI, 'filePathRaw', filePathRaw);
    set(handles.RawDataLabel, 'String', filePathRaw);
    setappdata(handles.GUI, 'filePathProc', filePathProc);
    set(handles.ProcDataLabel, 'String', filePathProc);
end



% --- Executes during object creation, after setting all properties.
function DateTextEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DateTextEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PipSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PipSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PipSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of PipSizeEdit as a double
pipSize = str2double(get(hObject,'String'))
setappdata(handles.GUI, 'pipSize', pipSize);


% --- Executes during object creation, after setting all properties.
function PipSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PipSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MeasPipBtn.
function MeasPipBtn_Callback(hObject, eventdata, handles)
% hObject    handle to MeasPipBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% display a frame in the video, tell user to click top and bottom of
% pipette

displayFigNum = getappdata(handles.GUI, 'displayFigNum');
fig = figure(displayFigNum);
clf;
frames = getappdata(handles.GUI, 'frames');
imshow(frames(:,:,1));

[x, y] = ginput(2);
pipSize = round(abs(y(2) - y(1)));
setappdata(handles.GUI, 'pipSize', pipSize);
set(handles.PipSizeEdit, 'String', num2str(pipSize));



function EmbryoNameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to EmbryoNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EmbryoNameEdit as text
%        str2double(get(hObject,'String')) returns contents of EmbryoNameEdit as a double
embryoNum = get(hObject, 'String');
setappdata(handles.GUI, 'embryoNum', embryoNum);


% --- Executes during object creation, after setting all properties.
function EmbryoNameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EmbryoNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in ManualMeasure.
function ManualMeasure_Callback(hObject, eventdata, handles)
% hObject    handle to ManualMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ManualMeasure
if get(hObject,'Value')
    manualMeasure = 1;
else
    manualMeasure = 0;
end

setappdata(handles.GUI, 'manualMeasure', manualMeasure);



% --- Executes on button press in EmbryoMeasBtn.
function EmbryoMeasBtn_Callback(hObject, eventdata, handles)
% hObject    handle to EmbryoMeasBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frames = getappdata(handles.GUI, 'frames');
t = getappdata(handles.GUI, 't');
params = getappdata(handles.GUI, 'params');
embryoNum = getappdata(handles.GUI, 'embryoNum');
manualMeasure = getappdata(handles.GUI, 'manualMeasure');
cannyThresh = params.cannyThresh;
filePathRaw = getappdata(handles.GUI, 'filePathRaw');
filePathProc = getappdata(handles.GUI, 'filePathProc');
procFileName = params.procFileName;
extraFigNum = getappdata(handles.GUI, 'extraFigNum');
startFrame = getappdata(handles.GUI, 'currFrame');
alreadyCropped = getappdata(handles.GUI, 'alreadyCropped');

axes(handles.MeasAxes);
fig = figure(extraFigNum);

% 1. Get ROI around just pipette opening
[ROIframes] = GetPipetteROI(frames, cannyThresh, fig, filePathRaw, alreadyCropped);

if exist([filePathRaw '\pipRef.mat'], 'file')
    set(handles.PipRefIndicator, 'String', 'YES');
    setappdata(handles.GUI, 'pipRefExists', 1);
end

axes(handles.PlotAxes);
cla;
sROI = size(ROIframes);
params.sROI = sROI;
setappdata(handles.GUI, 'params', params);
setappdata(handles.GUI, 'frames', ROIframes);
setappdata(handles.GUI, 'extraFigNum', extraFigNum);
setappdata(handles.GUI, 'alreadyCropped', 1);
clear frames;

figure(extraFigNum); clf;

% 2. Measure params from aspiration depth
[paramsFit, fig] = MeasureEmbryoAspiration(ROIframes, t, params, embryoNum, ...
    manualMeasure, filePathProc, procFileName, handles, fig, startFrame);

setappdata(handles.GUI, 'paramsFit', paramsFit);

% 3. Display in GUI
axes(handles.PlotAxes);
set(handles.k1_box, 'String', sprintf('%0.3f',paramsFit.k1));
set(handles.n1_box, 'String', sprintf('%0.3f',paramsFit.n1));
set(handles.tau_box, 'String', sprintf('%0.3f',paramsFit.tau));
set(handles.k0_box, 'String', sprintf('%0.3f',paramsFit.k0));


% --- Executes on button press in LoadVideoBtn.
function LoadVideoBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LoadVideoBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDate = getappdata(handles.GUI, 'currDate')
dateU = getappdata(handles.GUI, 'dateU')
measType = getappdata(handles.GUI, 'measType')
embryoNum = getappdata(handles.GUI, 'embryoNum')
pipSize = getappdata(handles.GUI, 'pipSize')
manualMeasure = getappdata(handles.GUI, 'manualMeasure')
filePathRaw = getappdata(handles.GUI, 'filePathRaw')
filePathProc = getappdata(handles.GUI, 'filePathProc')

if isequal(filePathRaw, '') || isequal(filePathProc, '')
    errordlg('Please select raw and processed data locations!');
    return;
end

[frames, t, params] = LoadVideo(measType, currDate, pipSize, filePathRaw, ...
    embryoNum, handles);

currFrame = round(params.frameStartMult*params.frameRate);
setappdata(handles.GUI, 'frames', frames);
setappdata(handles.GUI, 't', t);
setappdata(handles.GUI, 'params', params);
setappdata(handles.GUI, 'videoLoaded', 1);
setappdata(handles.GUI, 'lastFrame', size(frames,3));
setappdata(handles.GUI, 'currFrame', currFrame);
setappdata(handles.GUI, 'alreadyCropped', 0);

% adjust slider params
set(handles.FrameSlider, 'Min', 1);
set(handles.FrameSlider, 'Max', size(frames,3));
set(handles.FrameSlider, 'Value', currFrame);

% reset some variables
setappdata(handles.GUI, 'zpEnter', 0);
setappdata(handles.GUI, 'cellEnter', 0);
axes(handles.PlotAxes);
cla;
set(handles.ZPframeLabel, 'String', 'Frame: 0');
set(handles.CellFrameLabel, 'String', 'Frame: 0');

set(handles.MeasPipBtn, 'Enable', 'on');
set(handles.EmbryoMeasBtn, 'Enable', 'on');
set(handles.NewDisplayBtn, 'Enable', 'on');
set(handles.NewDisplayBtnPlot, 'Enable', 'off');
set(handles.currFrameLabel, 'String', ['frame: ' ...
    num2str(round(params.frameStartMult*params.frameRate))]);
set(handles.videoFileLabel, 'String', ['file: E' ...
    num2str(embryoNum) '.avi']);

% --- Executes on button press in RawDataBtn.
function RawDataBtn_Callback(hObject, eventdata, handles)
% hObject    handle to RawDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.TypeSelector, 'String'));
measType = contents{get(handles.TypeSelector,'Value')};
filePathRaw = uigetdir(['C:\Users\Livia\Desktop\IVF\Raw Data\Videos\' measType], ...
    'Select Raw Data Folder');
setappdata(handles.GUI, 'filePathRaw', filePathRaw);
set(handles.RawDataLabel, 'String', filePathRaw);

% check for existence of pipette reference file
% do this every time raw data is changed
if exist([filePathRaw '\pipRef.mat'], 'file')
    set(handles.PipRefIndicator, 'String', 'YES');
    setappdata(handles.GUI, 'pipRefExists', 1);
else
    set(handles.PipRefIndicator, 'String', 'NO');
    setappdata(handles.GUI, 'pipRefExists', 0);
end


% --- Executes on button press in ProcDataBtn.
function ProcDataBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ProcDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.TypeSelector, 'String'));
measType = contents{get(handles.TypeSelector,'Value')};
filePathProc = uigetdir(['C:\Users\Livia\Desktop\IVF\Processed Data\' measType], ...
    'Select Processed Data Folder');
setappdata(handles.GUI, 'filePathProc', filePathProc);
set(handles.ProcDataLabel, 'String', filePathProc);


% --- Executes during object creation, after setting all properties.
function GUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in NewDisplayBtn.
function NewDisplayBtn_Callback(hObject, eventdata, handles)
% hObject    handle to NewDisplayBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


frames = getappdata(handles.GUI, 'frames');
displayFigNum = getappdata(handles.GUI, 'displayFigNum');

if ~isempty(frames)
    fig = figure(displayFigNum); clf;
    imshow(frames(:,:,1));
else
    errordlg('Error: no frames loaded');
end



function k1_box_Callback(hObject, eventdata, handles)
% hObject    handle to k1_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k1_box as text
%        str2double(get(hObject,'String')) returns contents of k1_box as a double


% --- Executes during object creation, after setting all properties.
function k1_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k1_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n1_box_Callback(hObject, eventdata, handles)
% hObject    handle to n1_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n1_box as text
%        str2double(get(hObject,'String')) returns contents of n1_box as a double


% --- Executes during object creation, after setting all properties.
function n1_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n1_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_box_Callback(hObject, eventdata, handles)
% hObject    handle to tau_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau_box as text
%        str2double(get(hObject,'String')) returns contents of tau_box as a double


% --- Executes during object creation, after setting all properties.
function tau_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k0_box_Callback(hObject, eventdata, handles)
% hObject    handle to k0_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k0_box as text
%        str2double(get(hObject,'String')) returns contents of k0_box as a double


% --- Executes during object creation, after setting all properties.
function k0_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k0_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NewDisplayBtnPlot.
function NewDisplayBtnPlot_Callback(hObject, eventdata, handles)
% hObject    handle to NewDisplayBtnPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

paramsFit = getappdata(handles.GUI, 'paramsFit');
plotFigNum = getappdata(handles.GUI, 'plotFigNum');

if ~isempty(paramsFit)
    fig = figure(plotFigNum); clf;
    KelvinFit3(t, A, Fin, 1, [paramsFit.k0 paramsFit.k1 ...
        paramsFit.tau paramsFit.n1]);
else
    errordlg('Error: no params calculated');
end


% --- Executes on button press in RightBtn.
function RightBtn_Callback(hObject, eventdata, handles)
% hObject    handle to RightBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frames = getappdata(handles.GUI, 'frames');
currFrame = getappdata(handles.GUI, 'currFrame');
lastFrame = getappdata(handles.GUI, 'lastFrame');

if ~isempty(frames) && lastFrame > 0 && currFrame > 0
    if currFrame < lastFrame
        axes(handles.MeasAxes);
        imshow(frames(:,:,currFrame + 1));
        setappdata(handles.GUI, 'currFrame', currFrame + 1);
        set(handles.currFrameLabel, 'String', ['frame: ' ...
            num2str(currFrame + 1)]);
        set(handles.FrameSlider, 'Value', currFrame + 1);
    else
        errordlg('Already displaying last frame');
    end
else
    errordlg('Error: no frames loaded');
end

% --- Executes during object creation, after setting all properties.
function RightBtn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RightBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

rightArrow = imread('rightarrow.jpg');
set(hObject,'units','pixels');
buttonsize = get(hObject,'position');
imagesize = size(rightArrow);
newsize = ceil(imagesize(1:2)/(1.5*max(imagesize(1:2)./buttonsize([4 3]))));
newimage = rightArrow(round(linspace(1,imagesize(1),newsize(1))),...
    round(linspace(1,imagesize(2),newsize(2))),:);
set(hObject,'cdata',newimage);

% --- Executes on button press in LeftBtn.
function LeftBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LeftBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frames = getappdata(handles.GUI, 'frames');
currFrame = getappdata(handles.GUI, 'currFrame');
lastFrame = getappdata(handles.GUI, 'lastFrame');

if ~isempty(frames) && lastFrame > 0 && currFrame > 0
    if currFrame > 1
        axes(handles.MeasAxes);
        imshow(frames(:,:,currFrame - 1));
        setappdata(handles.GUI, 'currFrame', currFrame - 1);
        set(handles.currFrameLabel, 'String', ['frame: ' ...
            num2str(currFrame - 1)]);
        set(handles.FrameSlider, 'Value', currFrame - 1);
    else
        errordlg('Already displaying first frame');
    end
else
    errordlg('Error: no frames loaded');
end

% --- Executes during object creation, after setting all properties.
function LeftBtn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LeftBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

leftArrow = imread('leftarrow.jpg');
set(hObject,'units','pixels');
buttonsize = get(hObject,'position');
imagesize = size(leftArrow);
newsize = ceil(imagesize(1:2)/(1.5*max(imagesize(1:2)./buttonsize([4 3]))));
newimage = leftArrow(round(linspace(1,imagesize(1),newsize(1))),...
    round(linspace(1,imagesize(2),newsize(2))),:);
set(hObject,'cdata',newimage);


% % --- Executes on button press in PipRefBtn.
% function PipRefBtn_Callback(hObject, eventdata, handles)
% % hObject    handle to PipRefBtn (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% filePathRaw = getappdata(handles.GUI, 'filePathRaw');
% extraFig = getappdata(handles.GUI, 'extraFig');
% MakeRoiTemplate(filePathRaw, .35, extraFig);


% --- Executes on button press in ZonaEnterBtn.
function ZonaEnterBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ZonaEnterBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currFrame = getappdata(handles.GUI, 'currFrame');
setappdata(handles.GUI, 'zpEnter', currFrame);
set(handles.ZPframeLabel, 'String', ['Frame: ' num2str(currFrame)]);

% --- Executes on button press in CellEnterBtn.
function CellEnterBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CellEnterBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currFrame = getappdata(handles.GUI, 'currFrame');
setappdata(handles.GUI, 'cellEnter', currFrame);
set(handles.CellFrameLabel, 'String', ['Frame: ' num2str(currFrame)]);

% --- Executes on button press in SaveFrameBtn.
function SaveFrameBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFrameBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

zpEnter = getappdata(handles.GUI, 'zpEnter');
cellEnter = getappdata(handles.GUI, 'cellEnter');
filePathProc = getappdata(handles.GUI, 'filePathProc');
embryoNum = getappdata(handles.GUI, 'embryoNum');

if zpEnter > 0
    if ~isequal(filePathProc, '')
        save([filePathProc '\celldata_E' embryoNum '.mat'], 'zpEnter', ...
            'cellEnter');
    else
        errordlg('Error: no processed data path selected');
    end
else
    errordlg('Error: must select at least ZP frame');
end


% --- Executes on button press in LoadDataBtn.
function LoadDataBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filePathProc = getappdata(handles.GUI, 'filePathProc');
embryoNum = getappdata(handles.GUI, 'embryoNum');
dateU = getappdata(handles.GUI, 'dateU');
measType = getappdata(handles.GUI, 'measType');

if isequal(measType, 'Human')
    extraString = '_human';
else
    extraString = '';
end

testPath = [filePathProc '\AutoMeasure\aspiration_data_' dateU ...
    extraString '_E' embryoNum '.mat'];

if ~exist(testPath, 'file')
    testPath = [filePathProc '\aspiration_data_' dateU extraString ...
        '_E' embryoNum '.mat'];
end

if exist(testPath, 'file')
   
    load(testPath);
    axes(handles.PlotAxes);
    cla;
    
    L = length(aspiration_depth);
    plot(t, 10^6*[F0/(k0 + k1) A(1:(L-1))],'ob', 'Color', [0 0 1]);
    hold on;
    grid on;
    % this needs to be shifted to account for truncation during fitting
    plot(xfine, 10^6*yfit, 'Color', [0 0 1]);
    set(gca, 'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
    title('Aspiration Depth into Micropipette');
    
    % 3. Display in GUI
    set(handles.k1_box, 'String', sprintf('%0.3f',k1));
    set(handles.n1_box, 'String', sprintf('%0.3f',n1));
    set(handles.tau_box, 'String', sprintf('%0.3f',tau));
    set(handles.k0_box, 'String', sprintf('%0.3f',k0));
    
else
    errordlg('This embryo has not been measured yet!');
end


% --- Executes on slider movement.
function FrameSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FrameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

frames = getappdata(handles.GUI, 'frames');
value = get(hObject,'Value');
currFrame = getappdata(handles.GUI, 'currFrame');
lastFrame = getappdata(handles.GUI, 'lastFrame');

if ~isempty(frames)
        axes(handles.MeasAxes);
        imshow(frames(:,:, round(value)));
        setappdata(handles.GUI, 'currFrame', round(value));
        set(handles.currFrameLabel, 'String', ['frame: ' ...
            num2str(uint8(round(value)))]);
else
    errordlg('Error: no frames loaded');
end



% --- Executes during object creation, after setting all properties.
function FrameSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


