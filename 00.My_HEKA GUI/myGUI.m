function varargout = myGUI(varargin)
% MYGUI MATLAB code for myGUI.fig
%      MYGUI, by itself, creates a new MYGUI or raises the existing
%      singleton*.
%
%      H = MYGUI returns the handle to a new MYGUI or the handle to
%      the existing singleton*.
%
%      MYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYGUI.M with the given input arguments.
%
%      MYGUI('Property','Value',...) creates a new MYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before myGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to myGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help myGUI

% Last Modified by GUIDE v2.5 12-Jan-2021 20:38:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @myGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @myGUI_OutputFcn, ...
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


% --- Executes just before myGUI is made visible.
function myGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to myGUI (see VARARGIN)

% Choose default command line output for myGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes myGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = myGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile();

if ~(isnumeric(filename)&&isnumeric(pathname))
    
    filestring = strcat(pathname,filename);
    Data = HEKA_Importer(filestring);
    DataTable = Data.RecTable;
    handles.output.UserData = DataTable;
    
end

 RowNum = sum(table2array(handles.output.UserData(:,6))) + size(handles.output.UserData,1);
 SweepNum = table2array(handles.output.UserData(:,6));
 RecNum = table2array(handles.output.UserData(:,3));
 listStrings = cell(RowNum,1);
 
 indx = 0;
 for ii = 1:size(RecNum,1)
     indx = indx + 1;
     listStrings{indx} = strcat('Rec',num2str(RecNum(ii)));
     for jj = 1:SweepNum(ii)
         indx = indx + 1;
         listStrings{indx} = strcat('            Sweep ',num2str(jj));
     end
 end
 set(handles.listbox1,'String',listStrings);
 set(handles.text1,'String',filename);


% --------------------------------------------------------------------
function View_Callback(hObject, eventdata, handles)
% hObject    handle to View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 SweepNum = table2array(handles.output.UserData(:,6));
 RecNum = table2array(handles.output.UserData(:,3));RecIndx = zeros(length(RecNum),1);
 ListValue = get(handles.listbox1,'Value');
 for ii = 1:length(RecNum)
      if ii == 1
         RecIndx(ii) = 1;
      else
         RecIndx(ii) = RecIndx(ii-1) + SweepNum(ii-1)+1;
      end
 end
 
 if isempty(find(RecIndx == ListValue,1))

     RecSelect = length(RecIndx) - sum(ListValue < RecIndx); 
     SweepSelect = ListValue - RecIndx(RecSelect);
     RecUnit = table2array(handles.output.UserData(RecSelect,23));
     RecUnit = char(RecUnit{1});
     
     DataRaw = table2array(handles.output.UserData(RecSelect,21));
     DataRaw = cell2mat(DataRaw{1,1});
     DataStim = table2array(handles.output.UserData(RecSelect,22)); 
     DataStim = DataStim {1}.DA_3;
     if size(DataRaw,2) == size(DataStim,2)
        RawX = DataRaw(:,SweepSelect)*1000;
        StimX = DataStim(:,SweepSelect);
    end
    if (size(DataStim,2) == 1)&& (size(DataRaw,2) ~= size(DataStim,2))
        RawX = DataRaw(:,SweepSelect)*1000;
        StimX = DataStim(:,1);
    end
 
    f = table2array(handles.output.UserData(RecSelect,19)); % sampleing rate 40KHz
    time = (1:length(RawX))./f;
 
     axes(handles.axes1)
     cla(handles.axes1,'reset');zoom off;pan off;
     line(time,RawX);
     axis([0 time(end) floor((min(RawX)-20)/20)*20 ceil((max(RawX)+20)/20)*20])
     title('Raw Data');xlabel('Time (s)');
     if RecUnit == 'A'
         ylabel('Sweep (mV)');grid on;
     else
         ylabel('Sweep (pA)');grid on;
     end
     handles.axes1.UserData = [time',RawX];
     
     axes(handles.axes2)
     cla(handles.axes2,'reset');
     line(time,StimX);
     axis([0 time(end) round((min(StimX)-0.1),1) round((max(StimX)+0.1),1)])
     title('Stimulus Wave');xlabel('Time (s)');
     if RecUnit == 'A'
         ylabel('Stimulus Signal (pA)');grid on;
     else
         ylabel('Stimulus Signal (mV)');grid on;
     end
      handles.axes2.UserData = [time',StimX];
     
 else
     axes(handles.axes1);
     cla(handles.axes1,'reset');
     axes(handles.axes2);
     cla(handles.axes2,'reset');
     
 end
 


% --------------------------------------------------------------------
function Export_Callback(hObject, eventdata, handles)
% hObject    handle to Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 SweepNum = table2array(handles.output.UserData(:,6));
 RecNum = table2array(handles.output.UserData(:,3));RecIndx = zeros(length(RecNum),1);
 ListValue = get(handles.listbox1,'Value');
 for ii = 1:length(RecNum)
      if ii == 1
         RecIndx(ii) = 1;
      else
         RecIndx(ii) = RecIndx(ii-1) + SweepNum(ii-1)+1;
      end
 end
 
 if isempty(find(RecIndx == ListValue,1))

     RecSelect = length(RecIndx) - sum(ListValue < RecIndx); 
     SweepSelect = ListValue - RecIndx(RecSelect);
     RecUnit = table2array(handles.output.UserData(RecSelect,23));
     RecUnit = char(RecUnit{1});
     
     DataRaw = table2array(handles.output.UserData(RecSelect,21));
     DataRaw = cell2mat(DataRaw{1,1});
     DataStim = table2array(handles.output.UserData(RecSelect,22)); 
     DataStim = DataStim {1}.DA_3;
     if size(DataRaw,2) == size(DataStim,2)
        RawX = DataRaw(:,SweepSelect)*1000;
        StimX = DataStim(:,SweepSelect);
    end
    if (size(DataStim,2) == 1)&& (size(DataRaw,2) ~= size(DataStim,2))
        RawX = DataRaw(:,SweepSelect)*1000;
        StimX = DataStim(:,1);
    end
  
    save_path = uigetdir('Save to ',pwd);
    
    if ischar(save_path)
        filename = strcat('RawData_','Rec',num2str(RecSelect),'_Sweep',num2str(SweepSelect),'.txt');
        dlmwrite(strcat(save_path,'\',filename),RawX);
        hm = msgbox(strcat(filename,' was saved in ',{' '},save_path));
        
        pause(3);delete(hm);
    end
    
 end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Tool_Callback(hObject, eventdata, handles)
% hObject    handle to Tool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ZoomTool_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (~isempty(handles.axes1.UserData))||(~isempty(handles.axes2.UserData))
    pan off; 
    zoom on;
end

% --------------------------------------------------------------------
function PanTool_Callback(hObject, eventdata, handles)
% hObject    handle to PanTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (~isempty(handles.axes1.UserData))||(~isempty(handles.axes2.UserData))
      zoom off; 
      pan on;
end
