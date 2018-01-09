function varargout = FEA_Interface(varargin)
% FEA_INTERFACE MATLAB code for FEA_Interface.fig
%      FEA_INTERFACE, by itself, creates a new FEA_INTERFACE or raises the existing
%      singleton*.
%
%      H = FEA_INTERFACE returns the handle to a new FEA_INTERFACE or the handle to
%      the existing singleton*.
%
%      FEA_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FEA_INTERFACE.M with the given input arguments.
%
%      FEA_INTERFACE('Property','Value',...) creates a new FEA_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FEA_Interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FEA_Interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FEA_Interface

% Last Modified by GUIDE v2.5 11-Mar-2013 10:07:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FEA_Interface_OpeningFcn, ...
                   'gui_OutputFcn',  @FEA_Interface_OutputFcn, ...
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


% --- Executes just before FEA_Interface is made visible.
function FEA_Interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FEA_Interface (see VARARGIN)

% Choose default command line output for FEA_Interface

handles.nelx=100;
handles.nely=15;
handles.nelz=15;
handles.aspectab=10;
handles.aspectac=10;
handles.styecast_ratio=2;
handles.bcmode=1;
handles.ex=0;
handles.ey=6.62e-4;
handles.ez=1.83e-3;
handles.elevation=0.5;
handles.line_pos=0.5;
handles.strain_to_plot=1;
handles.hold='Min';
handles.plane=1;
handles.line=1;
handles.strain_str='\epsilon_{xx}';
handles.nodex=101;
handles.nodey=16;
handles.nodez=16;

handles.output = hObject;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FEA_Interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FEA_Interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function nelx_Callback(hObject, eventdata, handles)
% hObject    handle to nelx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nelx as text
%        str2double(get(hObject,'String')) returns contents of nelx as a double

handles.nelx=str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nelx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nelx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nely_Callback(hObject, eventdata, handles)
% hObject    handle to nely (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nely as text
%        str2double(get(hObject,'String')) returns contents of nely as a double

handles.nely=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nely_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nely (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nelz_Callback(hObject, eventdata, handles)
% hObject    handle to nelz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nelz as text
%        str2double(get(hObject,'String')) returns contents of nelz as a double

handles.nelz=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nelz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nelz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bcmode.
function bcmode_Callback(hObject, eventdata, handles)
% hObject    handle to bcmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns bcmode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from bcmode

val = get(hObject,'Value');
str = get(hObject, 'String');

switch str{val}
    
    case 'Only edges constrained'
        
        
    handles.bcmode=1;
    
    case 'Ends in styecast'
        
        handles.bcmode=2;
        
    case 'Only bottom constrained'
        
        handles.bcmode=3;
        
end
       
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function bcmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bcmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ex_Callback(hObject, eventdata, handles)
% hObject    handle to ex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ex as text
%        str2double(get(hObject,'String')) returns contents of ex as a double


handles.ex=str2double(get(hObject,'String'))/100;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function aspectab_Callback(hObject, eventdata, handles)
% hObject    handle to aspectab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aspectab as text
%        str2double(get(hObject,'String')) returns contents of aspectab as a double

handles.aspectab=str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function aspectab_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aspectab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function aspectac_Callback(hObject, eventdata, handles)
% hObject    handle to aspectac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aspectac as text
%        str2double(get(hObject,'String')) returns contents of aspectac as a double

handles.aspectac=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function aspectac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aspectac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ey_Callback(hObject, eventdata, handles)
% hObject    handle to ey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ey as text
%        str2double(get(hObject,'String')) returns contents of ey as a double

handles.ey=str2double(get(hObject,'String'))/100;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function ey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ez_Callback(hObject, eventdata, handles)
% hObject    handle to ez (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ez as text
%        str2double(get(hObject,'String')) returns contents of ez as a double

handles.ez=str2double(get(hObject,'String'))/100;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ez_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ez (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elevation_Callback(hObject, eventdata, handles)
% hObject    handle to elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elevation as text
%        str2double(get(hObject,'String')) returns contents of elevation as a double

handles.elevation=str2double(get(hObject,'String'))/100;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function elevation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sy_Callback(hObject, eventdata, handles)
% hObject    handle to sy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sy as text
%        str2double(get(hObject,'String')) returns contents of sy as a double

handles.sy=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sz_Callback(hObject, eventdata, handles)
% hObject    handle to sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sz as text
%        str2double(get(hObject,'String')) returns contents of sz as a double

handles.sz=str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in surface.
function surface_Callback(hObject, eventdata, handles)
% hObject    handle to surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% send user-inputted parameters for validation

validation_script;

file_handling_script;


[sx,sy,sz]=Outside_Surface(SX,SY,SZ);

% Make these into a Slice plot
 
 figure


slice(SX,SY,SZ,SSTRAIN(:,:,:,handles.strain_to_plot),sx,sy,sz);
shading interp;
axis equal;
axis vis3d;
colormap jet;
%colorbar;
xlabel('x/mm');
ylabel('y/mm');
zlabel('z/mm');
    
clear SX SY SZ SSTRAIN X Y Z Disp dimensions filename filepath    





% --- Executes on button press in line.
function line_Callback(hObject, eventdata, handles)
% hObject    handle to line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



validation_script;

file_handling_script;

 
if handles.line==1 && handles.plane==2
    
    error('Line direction chosen is not in plane selected');
    
elseif handles.line==2 && handles.plane==3
    
    error('Line direction chosen is not in plane selected');
    
elseif handles.line==3 &&handles.plane==1

     error('Line direction chosen is not in plane selected');


end


            minx=min(SX(1,:,1));
            maxx=max(SX(1,:,1));
            miny=min(SY(:,1,1));
            maxy=max(SY(:,1,1));
            minz=min(SZ(1,1,:));
            maxz=max(SZ(1,1,:));
            nodex=size(SX(1,:,1),2);
            nodey=size(SY(:,1,1),1);
            nodez=size(SZ(1,1,:),3);
switch handles.plane
    
    case 1 
        
        z_pos=minz+(maxz-minz)*handles.elevation;
        
        if handles.line==1
            
            
           no_points=nodex;
           y_pos=miny+(maxy-miny)*handles.line_pos;
           
            
            
           xi=linspace(minx,maxx,no_points);
           yi=repmat(y_pos,1,no_points);
           zi=repmat(z_pos,1,no_points);
           
           VI=interp3(SX,SY,SZ,SSTRAIN(:,:,:,handles.strain_to_plot),xi,yi,zi);
figure
plot(xi,VI), axis square, xlabel('x'), ylabel( handles.strain_str);
           
        elseif handles.line==2
            
           no_points=nodey;
           x_pos=minx+(maxx-minx)*handles.line_pos;
           
           xi=repmat(x_pos,1,no_points);
           yi=linspace(miny,maxy,no_points);
           zi=repmat(z_pos,1,no_points);
           
           VI=interp3(SX,SY,SZ,SSTRAIN(:,:,:,handles.strain_to_plot),xi,yi,zi);
figure
plot(yi,VI), axis square, xlabel('y'), ylabel( handles.strain_str);
           
        end
        
    case 2 
        
        x_pos=minx+(maxx-minx)*handles.elevation;
        
        if handles.line==2
            
            no_points=handles.nodey;
            z_pos=minz+(maxz-minz)*handles.elevation;
            
            xi=repmat(x_pos,1,no_points);
            yi=linspace(miny,maxy,no_points);
            zi=repmat(z_pos,1,no_points);
            
            VI=interp3(SX,SY,SZ,SSTRAIN(:,:,:,handles.strain_to_plot),xi,yi,zi);
figure
plot(yi,VI), axis square, xlabel('y'), ylabel( handles.strain_str);
           
        elseif handles.line==3
            
            no_points=handles.nodez;
            y_pos=miny+(maxy-miny)*handles.line_pos;
            
            xi=repmat(x_pos,1,no_points);
            yi=repmat(y_pos,1,no_points);
            zi=linspace(minz,maxz,no_points);
            
            VI=interp3(SX,SY,SZ,SSTRAIN(:,:,:,handles.strain_to_plot),xi,yi,zi);
figure
plot(zi,VI), axis square, xlabel('z'), ylabel( handles.strain_str);
            
        end
        
    case 3 
        
        y_pos=miny+(maxy-miny)*handles.elevation;
        
        if handles.line==1
            
            no_points=handles.nodex;
            z_pos=minz+(maxz-minz)*handles.line_pos;
            
            xi=linspace(minx,maxx,no_points);
            yi=repmat(y_pos,1,no_points);
            zi=repmat(z_pos,1,no_points);
            
            VI=interp3(SX,SY,SZ,SSTRAIN(:,:,:,handles.strain_to_plot),xi,yi,zi);
figure
plot(xi,VI), axis square, xlabel('x'), ylabel( handles.strain_str);
            
        elseif handles.line==3
            
            no_points=handles.nodez;
            x_pos=minx+(maxx-minx)*handles.line_pos;
            
            xi=repmat(x_pos,1,no_points);
            yi=repmat(y_pos,1,no_points);
            zi=linspace(minz,maxz,no_points);
            
            VI=interp3(SX,SY,SZ,SSTRAIN(:,:,:,handles.strain_to_plot),xi,yi,zi);
figure
plot(zi,VI), axis square, xlabel('z'), ylabel( handles.strain_str);
        end
        
end



        
      
            
        
     

% --- Executes on button press in slice.
function slice_Callback(hObject, eventdata, handles)
% hObject    handle to slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


validation_script;

file_handling_script;

switch handles.plane
    
    case 1
elev_val=min(SZ(1,1,:))+(max(SZ(1,1,:))-min(SZ(1,1,:)))*handles.elevation;

sx=inf;
sy=inf;
sz=elev_val;

    case 2
        
        
        elev_val=min(SX(1,:,1))+(max(SX(1,:,1))-min(SX(1,:,1)))*handles.elevation;

sx=elev_val;
sy=inf;
sz=inf;


    case 3
        
        elev_val=min(SY(:,1,1))+(max(SY(:,1,1))-min(SY(:,1,1)))*handles.elevation;

sx=inf;
sy=elev_val;
sz=inf;

end
% Make these into a Slice plot
 
 figure


slice(SX,SY,SZ,SSTRAIN(:,:,:,handles.strain_to_plot),sx,sy,sz);
shading interp;
axis equal;
axis vis3d;
colormap jet;
%colorbar;
xlabel('x/mm');
ylabel('y/mm');
zlabel('z/mm');
    
clear SX SY SZ SSTRAIN X Y Z Disp dimensions filename filepath   


function styecast_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to styecast_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of styecast_ratio as text
%        str2double(get(hObject,'String')) returns contents of styecast_ratio as a double


styecast_ratio=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function styecast_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to styecast_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hold.
function hold_Callback(hObject, eventdata, handles)
% hObject    handle to hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hold

handles.hold=get(hObject,'Value')
guidata(hObject, handles);


% --- Executes on selection change in plane_to_plot.
function plane_to_plot_Callback(hObject, eventdata, handles)
% hObject    handle to plane_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plane_to_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plane_to_plot

val = get(hObject,'Value');
str = get(hObject, 'String');

switch val
    
    case 1
        
        
    handles.plane=1;
    set(handles.elevation_label, 'String', 'z');
    
    case 2
        
        handles.plane=2;
        
set(handles.elevation_label, 'String', 'x');
        
    case 3
        
        handles.plane=3;
        set(handles.elevation_label, 'String', 'y');
        
end
       
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function plane_to_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plane_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function line_pos_Callback(hObject, eventdata, handles)
% hObject    handle to line_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of line_pos as text
%        str2double(get(hObject,'String')) returns contents of line_pos as a double
handles.line_pos=str2double(get(hObject,'String'))/100;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function line_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to line_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in line_to_plot.
function line_to_plot_Callback(hObject, eventdata, handles)
% hObject    handle to line_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns line_to_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from line_to_plot

val = get(hObject,'Value');
str = get(hObject, 'String');

switch str{val}
    
    case 'x'
        
   handles.line=1;
   
   switch handles.plane
       
       case 1
           set(handles.line_label, 'String', 'y');
           
       case 2 
          error('Line direction does not lie in plane selected');
          
       case 3
           set(handles.line_label, 'String', 'z');
           
           
   end
    case 'y'
        
        handles.line=2;
        
           switch handles.plane
       
       case 1
           set(handles.line_label, 'String', 'x');
           
       
       case 2
           set(handles.line_label, 'String', 'z');
           
       case 3
          error('Line direction does not lie in plane selected');
           end

        
    case 'z'
        
        handles.line=3;
        
        switch handles.plane
        
        case 1
          error('Line direction does not lie in plane selected');
       
        case 2
           set(handles.line_label, 'String', 'y');
           
        
        case 3
           set(handles.line_label, 'String', 'x');
           
        end
       
       
        
        
        
end
       
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function line_to_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to line_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in strain_to_plot.
function strain_to_plot_Callback(hObject, eventdata, handles)
% hObject    handle to strain_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns strain_to_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from strain_to_plot

val = get(hObject,'Value');
str = get(hObject, 'String');

switch str{val}
    
    case 'xx'
        
        
    handles.strain_to_plot=1;
    handles.strain_str='\epsilon_{xx}';
    
    case 'yy'
        
        handles.strain_to_plot=2;
        handles.strain_str='\epsilon_{yy}';

        
    case 'zz'
        
        handles.strain_to_plot=3;
        handles.strain_str='\epsilon_{zz}';
    case 'xy'
        
        handles.strain_to_plot=4;
        handles.strain_str='\epsilon_{xy}';
    case 'zy'
        
        handles.strain_to_plot=5;
        handles.strain_str='\epsilon_{zy}';
    case 'zx'
        
        handles.strain_to_plot=6;
        handles.strain_str='\epsilon_{zx}';
end
       
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function strain_to_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strain_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
