function varargout = EmpPatterns(varargin)
% EmpPatterns M-file for EmpPatterns.fig
%   EmpPatterns, by itself, creates a new EmpPatterns or raises the existing
%   singleton*.
%
%   H = EmpPatterns returns the handle to a new EmpPatterns or the handle to
%   the existing singleton*.
%
%   EmpPatterns('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in EmpPatterns.m with the given input arguments.
%
%   EmpPatterns('Property','Value',...) creates a new EmpPatterns or raises the
%   existing singleton*. Starting from the left, property value pairs are
%   applied to the GUI before EmpPatterns_OpeningFcn gets called. An
%   unrecognized property name or invalid value makes property application
%   stop. All inputs are passed to EmpPatterns_OpeningFcn via varargin.
%
%   *See GUI Options on GUIDE's Tools menu. Choose "GUI allows only one
%   instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton', gui_Singleton, ...
    'gui_OpeningFcn', @EmpPatterns_OpeningFcn, ...
    'gui_OutputFcn',  @EmpPatterns_OutputFcn, ...
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


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes just before EmpPatterns is made visible.
function EmpPatterns_OpeningFcn2(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EmpPatterns (see VARARGIN)

% Choose default command line output for EmpPatterns
handles.output = hObject;
% Initialize handles for UI components
handles.PatternStyle = findobj('Tag', 'PatternStyle');
handles.PatternScheme = findobj('Tag', 'PatternScheme');
handles.n1 = findobj('Tag', 'n1');
handles.k1 = findobj('Tag', 'k1');
handles.r1 = findobj('Tag', 'r1');
handles.ch1 = findobj('Tag', 'ch1');
handles.text5 = findobj('Tag', 'text5');
% Update handles structure
guidata(hObject, handles);

% --- Executes just before EmpPatterns is made visible.
function EmpPatterns_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EmpPatterns (see VARARGIN)

% Choose default command line output for EmpPatterns
handles.output = hObject;
handles.samples = [];  % Add this line to initialize samples field
handles.datasize = [];  % Add this line to initialize datasize field
handles.trialcount = [];  % Add this line to initialize trialcount field
handles.TableTitle = [];
% Update handles structure
guidata(hObject, handles);

% Wait for the GUI to be fully constructed before attempting to access handles
% uiwait(hObject);

% Now, retrieve handles using findobj
handles = guihandles(hObject);

% Save the updated handles structure
guidata(hObject, handles);




% UIWAIT makes EmpPatterns wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = EmpPatterns_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;

% --- Executes on selection change in PatternStyle.
function PatternStyle_Callback(hObject, eventdata, handles)
% hObject    handle to PatternStyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns PatternStyle contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PatternStyle

% --- Executes during object creation, after setting all properties.
function PatternStyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PatternStyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in PatternScheme.
function PatternScheme_Callback(hObject, eventdata, handles)
% hObject    handle to PatternScheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns PatternScheme contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PatternScheme

% --- Executes during object creation, after setting all properties.
function PatternScheme_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PatternScheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k1_Callback(hObject, eventdata, handles)
% hObject    handle to k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k1 as text
%        str2double(get(hObject,'String')) returns contents of k1 as a double

% --- Executes during object creation, after setting all properties.
function k1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r1_Callback(hObject, eventdata, handles)
% hObject    handle to r1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r1 as text
%        str2double(get(hObject,'String')) returns contents of r1 as a double

% --- Executes during object creation, after setting all properties.
function r1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n1_Callback(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n1 as text
%        str2double(get(hObject,'String')) returns contents of n1 as a double


% --- Executes during object creation, after setting all properties.
function n1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch1_Callback(hObject, eventdata, handles)
% hObject    handle to ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1 as text
%        str2double(get(hObject,'String')) returns contents of ch1 as a double

% --- Executes during object creation, after setting all properties.
function ch1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global samples;
samples = load(uigetfile('*.csv', 'Select the data file containing'));
[nrows,ncols]=size(samples);
set(handles.datasize, 'String', nrows);
set(handles.trialcount, 'String', ncols);
set(handles.n1, 'String', ncols);



% sequence'));
% Ask the user to select a text file


% global samples;
% [filename, pathname] = uigetfile('*.txt', 'Select Text File');
% if isequal(filename, 0) || isequal(pathname, 0)
%     % User canceled the operation
%     return;
% end
%
% try
%     % Assuming the file contains a matrix of sample sequences
%     samples = load(fullfile(pathname, filename));
%
%     % You can update your GUI based on the loaded data if needed
%     % For example, update a listbox, display the matrix, etc.
%
%     % Update handles structure
%     handles.loadedData = samples;
%     guidata(hObject, handles);
%     samples=readtable(uigetfile('*.txt', 'Select the data file containing sequence'));
%
%     % Notify the user that the data has been successfully loaded
%     msgbox('Data loaded successfully!', 'Success', 'modal');
%
% catch exception
%     % Handle any errors that might occur during file reading
%     errordlg(['Error loading data from the file: ' exception.message], 'Error', 'modal');
% end


% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global n;global ch;
global samples;
% samples=load(uigetfile('*.csv', 'Select the data file containing sequence'));
val1=get(handles.PatternStyle,'value');
val2=get(handles.PatternScheme,'value');
n=str2num(get(handles.n1,'String'));
k=str2num(get(handles.k1,'String'));
r=str2num(get(handles.r1,'String'));
ch=str2num(get(handles.ch1,'String'));
c=0;
if length(k)*length(n)==0
    msgbox('  Either you have not entered n or k                 ', 'Error: Parameter values are not entered       ');
elseif val2==3&&length(k)*length(r)==0 || val2==4&&length(k)*length(r)==0
    msgbox('  Either you have not entered r or k                 ', 'Error: Parameter values are not entered       ');
else
    switch val1
        case 1
            switch val2
                case 1
                    table=Patterns_AMN(samples,ch,n,k,c);%Patterns_AMN(samples, ch, n, k, c)
                case 2
                    table=Patterns_AMM(samples,ch,n,k,c);
                case 3
                    table=Patterns_AMT(samples,ch,k,r,c);
                case 4
                    table=Patterns_AMW(samples,ch,k,r,c);
            end
        case 2
            switch val2
                case 1
                    table=Patterns_EN(samples,ch,n,k,c);
                case 2
                    table=Patterns_EM(samples,ch,n,k,c);
                case 3
                    table=Patterns_ET(samples,ch,k,r,c);
                case 4
                    table=Patterns_EW(samples,ch,k,r,c);
            end
        case 3
            switch val2
                case 1
                    table=Patterns_ALN(samples,ch,n,k,c);
                case 2
                    table=Patterns_ALM(samples,ch,n,k,c);
                case 3
                    table=Patterns_ALT(samples,ch,k,r,c);
                case 4
                    table=Patterns_ALW(samples,ch,k,r,c);

            end
        case 4
            switch val2
                case 1
                    table=Patterns_UN(samples,n,k,r);
                case 2
                    table=Patterns_UM(samples,n,k,r);
                case 3
                    table=Patterns_UT(samples,n,k,r);
                case 4
                    table=Patterns_UW(samples,n,k,r);
            end
    end
    disp(table);
    dlmwrite('output_table.txt', table, 'delimiter', '\t');

    if length(k) == 1
        table = tabulate(table);
        % Add the CDF to the data table
        table(:, 3) = cumsum(table(:, 2))./sum(table(:, 2));
        % Assuming you have a table called 'table' containing your results
        data = num2cell([table(:, 1:2), cumsum(table(:, 2))./sum(table(:, 2))]);
        % Get the position of "text19"
        texttitle_position = get(handles.text5, 'Position');
        column_names = {'X', 'Frequency', 'Empirical CDF'};
        % Create uitable
        uitable_handle = uitable('Parent', handles.figure1, 'Data', data, ...
            'Position', [texttitle_position(1) + 40, texttitle_position(2) + 70, 320, 200]);

        % Update the plot
        newFigure = figure;
        plot(table(:, 1), table(:, 3));
        ylim([0, 1]);
        xlim([min(table(:, 1)), max(table(:, 1))]);
        title('Empirical Probability Distribution');
    else
        texttitle_position = get(handles.text5, 'Position');
        set(handles.text5,'String','Correlation Matrix');
        uitable_handle = uitable('Parent', handles.figure1, 'Data', corrcoef(table), ...
            'Position', [texttitle_position(1) + 40, texttitle_position(2) + 70, 320, 200]);
        
        
        % Create a new figure for the plot matrix
        newFigure = figure;

        % Assuming 'table' contains your data
        plotmatrix(table);

        % Update the title
        title('Plot Matrix');
        
        
    end
end