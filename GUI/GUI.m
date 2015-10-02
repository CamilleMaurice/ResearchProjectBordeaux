function varargout = GUI(varargin)
%GUI M-file for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('Property','Value',...) creates a new GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI('CALLBACK') and GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 22-Sep-2015 16:03:30
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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
end

% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;


%Share data among callbacks 
data = struct('initialFrameName','default','OtherData',0);
set(handles.Openfile,'UserData',data);

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in varCheckBox.
function varCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to varCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of varCheckBox
end

% --- Executes on button press in maxflowCheckBox.
function maxflowCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to maxflowCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
% Hint: get(hObject,'Value') returns toggle state of maxflowCheckBox


% --------------------------------------------------------------------
function Openfile_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile({'*.png';'*.bmp';'*.jpeg'},'Load image');
% hObject    handle to Openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(hObject,'UserData');

fullImageFileName = fullfile(PathName, FileName);

data.initialFrameName = fullImageFileName;

% Store data in UserData of Openfile
set(hObject,'UserData',data);


im_original = imread(fullImageFileName);
image(im_original); 
hold on;
%pts = readPoints(image,2);
%rect = getrect(im_original);
%rect
h=imrect;
rect=uint16(getPosition(h));
rect;
end



function pts = readPoints(image, n)
%readPoints   Read manually-defined points from image
%   POINTS = READPOINTS(IMAGE) displays the image in the current figure,
%   then records the position of each click of button 1 of the mouse in the
%   figure, and stops when another button is clicked. The track of points
%   is drawn as it goes along. The result is a 2 x NPOINTS matrix; each
%   column is [X; Y] for one point.
% 
%   POINTS = READPOINTS(IMAGE, N) reads up to N points only.

if nargin < 2
    n = Inf;
    pts = zeros(2, 0);
else
    pts = zeros(2, n);
end

imshow(image);     % display image
xold = 0;
yold = 0;
k = 0;
hold on;           % and keep it there while we plot

while 1
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    k = k + 1;
    pts(1,k) = xi;
    pts(2,k) = yi;

      if xold
          plot([xold xi], [yold yi], 'go-');  % draw as we go
      else
          plot(xi, yi, 'go');         % first point on its own
      end

      if isequal(k, n)
          break
      end
      xold = xi;
      yold = yi;
  end

hold off;
if k < size(pts,2)
    pts = pts(:, 1:k);
end

end

% --- Executes on button press in StartReadSeq.
function StartReadSeq_Callback(hObject, eventdata, handles)
% hObject    handle to StartReadSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of StartReadSeq
% Get UserData from the Openfile
data = get(handles.Openfile,'UserData');
FileName = data.initialFrameName;
BaseFileName = FileName(1:end-5);
Ext = FileName(end-3:end);
%Read image sequence
for i=0:9       
    i
    sprintf('%s%d%s',BaseFileName,i,Ext);
    img = imread(sprintf('%s%d%s',BaseFileName,i,Ext) );    
    image(img);
    pause(1)
end
end