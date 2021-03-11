function varargout = project2(varargin)
% PROJECT2 MATLAB code for project2.fig
%      PROJECT2, by itself, creates a new PROJECT2 or raises the existing
%      singleton*.
%
%      H = PROJECT2 returns the handle to a new PROJECT2 or the handle to
%      the existing singleton*.
%
%      PROJECT2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECT2.M with the given input arguments.
%
%      PROJECT2('Property','Value',...) creates a new PROJECT2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before project2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to project2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help project2

% Last Modified by GUIDE v2.5 16-Jan-2021 01:47:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @project2_OpeningFcn, ...
    'gui_OutputFcn',  @project2_OutputFcn, ...
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


% --- Executes just before project2 is made visible.
function project2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to project2 (see VARARGIN)

% Choose default command line output for project2
handles.output = hObject;
fileID_opening = fopen('output.txt','w');
fclose(fileID_opening);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes project2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = project2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function equation_Callback(hObject, eventdata, handles)
% hObject    handle to equation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of equation as text
%        str2double(get(hObject,'String')) returns contents of equation as a double


% --- Executes during object creation, after setting all properties.
function equation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to equation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% % % % % % % % % % % % % % % % % table =uitable('columnname','equation','position',[20.8 400.769 300.2 170.615]);

% --- Executes on button press in insert.
function insert_Callback(hObject, eventdata, handles)
eq = get(handles.equation,'String');
row ={eq};
oldData = get(handles.eq_box,'String');
newData=[oldData; row];
set(handles.eq_box,'String',newData);
set(handles.equation,'String','');





% set(handles.equations,'String',get(handles.equation,'String'));

% hObject    handle to insert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in choose_methode.
function choose_methode_Callback(hObject, eventdata, handles)
% hObject    handle to choose_methode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns choose_methode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from choose_methode


% --- Executes during object creation, after setting all properties.
function choose_methode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choose_methode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
cla reset ;
set(handles.eq_box, 'String', '');
set(handles.equation, 'String', '');
set(handles.initial_box,'String','');
set(handles.max_iteration,'String','50');
set(handles.epsilon,'String','0.00001');
set(handles.choosen_file,'String','');
set(handles.init,'String','');
set(handles.answer_box1, 'String', '');
set(handles.answer_box2, 'String', '');
set(handles.answer_box3, 'String', '');
set(handles.answer_box4, 'String', '');
set(handles.time1, 'String', '');
set(handles.time2, 'String', '');
set(handles.time3, 'String', '');
set(handles.time4, 'String', '');
set(handles.iterations, 'String', '');
set(handles.precision, 'String', '');
set(handles.itable,'Data','');
delete output.txt ;

% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in choose_file.
function choose_file_Callback(hObject, eventdata, handles)
[filename pathname] = uigetfile({'*.txt'},'File Selector');
if(filename==0)
    set(handles.choosen_file, 'String', '');
else
    set(handles.eq_box, 'String', '');
    set(handles.initial_box, 'String', '');
    set(handles.choosen_file, 'String', strcat(pathname,filename));
    
    fid1=fopen(strcat(pathname,filename));
    numOfEquations = str2num(fgetl(fid1))
    methodName= fgetl(fid1)
    if strcmp(methodName,'Gaussian-elimination') || strcmp(methodName,'Gaussian-Elimination')
        set(handles.methods, 'Value', 1)
    end
    if strcmp(methodName,'LU decomposition')
        set(handles.methods, 'Value', 2)
    end
    if strcmp(methodName,'Gaussian-Jordan') || strcmp(methodName,'Gaussian-jordan')
        set(handles.methods, 'Value', 3)
    end
    if strcmp(methodName,'Gauss-Seidel') || strcmp(methodName,'Gauss-seidel')
        set(handles.methods, 'Value', 4)
    end
    
    
    
    for i=1:numOfEquations
        eq = fgetl(fid1);
        row ={eq };
        oldData = get(handles.eq_box,'String');
        newData=[oldData; row];
        set(handles.eq_box,'String',newData);
    end
    
    if( ~feof(fid1))
        init=fgetl(fid1);
        init = strsplit(init) %%may use str2double()
        set(handles.initial_box,'String',init);
    end
    fclose(fid1);
    
end
% hObject    handle to choose_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_iteration_Callback(hObject, eventdata, handles)
% hObject    handle to max_iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_iteration as text
%        str2double(get(hObject,'String')) returns contents of max_iteration as a double


% --- Executes during object creation, after setting all properties.
function max_iteration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epsilon_Callback(hObject, eventdata, handles)
% hObject    handle to epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilon as text
%        str2double(get(hObject,'String')) returns contents of epsilon as a double


% --- Executes during object creation, after setting all properties.
function epsilon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Gaussian_elimination.
function Gaussian_elimination_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussian_elimination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Gaussian_elimination


% --- Executes on button press in LU_decomposition.
function LU_decomposition_Callback(hObject, eventdata, handles)
% hObject    handle to LU_decomposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LU_decomposition


% --- Executes on button press in Gaussian_Jordan.
function Gaussian_Jordan_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussian_Jordan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Gaussian_Jordan


% --- Executes on button press in Gauss_Seidel.
function Gauss_Seidel_Callback(hObject, eventdata, handles)
% hObject    handle to Gauss_Seidel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Gauss_Seidel


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
selected_index=get(handles.methods,'Value')
Data = get(handles.eq_box,'String');
initials = get(handles.initial_box,'String'); 

if(isempty(get(handles.eq_box,'String')))
    errordlg('Please enter the equations','Error');
    return
end

len=length(Data);

syms x
syms y
syms z

syms o
syms m
syms n

syms h
syms f
syms g

char_arr=[];

try
    arr_final=[];
    result_arr = zeros(10,1);
    for k=1 : len
        d=Data(k);
        ss = [d{:}];
        [c,t] = coeffs(ss,[x y z o m])
       
        h=size(c);
        char_arr_size=size(char_arr);
        if char_arr_size(2) < h(2)
            char_arr = t 
        end
        arr=zeros(1,9);
        for i=1:h(2)
            if(t(i) == 'x')
                arr(1)=c(i);
            elseif(t(i) == 'y')
                arr(2)=c(i);
            elseif(t(i) == 'z')
                arr(3) =c(i);
            elseif(t(i) == 'o')
                arr(4) = c(i);
            elseif(t(i)== 'm')
                arr(5) = c(i);
            elseif(t(i)=='n')
                arr(6) = c(i);
            elseif(t(i) == 'h')
                arr(7) = c(i);
            elseif(t(i)== 'f')
                arr(8) = c(i);
            elseif(t(i) == 'g')
                arr(9) = c(i);
            elseif(t(i) == 1 )
                result_arr(k) = -1*c(i);
            end
        end
        arr_final=[arr_final;arr];
    end
    
    arr_size = size(arr_final);
    arr_size(1)
    
    result_arr
    arr_final
    
catch
    errordlg('May be incorect variable in equations','Error');
end

%%-----------------------------Gaussian-elimination-------------------
if(find(selected_index==1))
    now1=tic();
    n=len;
    A=arr_final(1:n,1:n);
    b=result_arr(1:n,:);
    x=zeros(n,1);
    for i = 1:n-1
        m= A(i+1:n,i)/A(i,i);
        if isnan(m)
            fprintf('ERROR no guasss elimination');
            errordlg('ERROR no guasss elimination','Error');
            return;
        elseif isinf(m)
            fprintf('ERROR no guasss elimination');
            errordlg('ERROR no guasss elimination','Error');
            return;
        else
            A(i+1:n,:) = A(i+1:n,:) - m*A(i,:);
            b(i+1:n,:) = b(i+1:n,:) - m*b(i,:);
        end
        
    end
    if A(end)==0 && b(end)~=0
        fprintf('solution doesnot exists');
        errordlg('solution doesnot exists','Error');
        return;
    end
    if b(end)==0
        fprintf('infinite number of solutios');
        errordlg('infinite number of solutios','Error');
        return;
    end
    
    x(n,:) = b(n,:)/A(n,n);
    for i = n-1:-1:1
        x(i,:) = (b(i,:) - A(i,i+1:n)*x(i+1:n,:))/A(i,i);
    end
    double(x);
    
    set(handles.answer_box1,'String',double(x));
    execution=toc(now1);
    set(handles.time1,'string',[num2str(execution) ' seconds'])
    
    fileID = fopen('output.txt','a');
    fprintf(fileID,'%s\n','Gaussian-elimination roots');
    fprintf(fileID,'%s     %s\n','var','value');
    for itr=1:len
        fprintf(fileID,'%s     %f\n',char(char_arr(itr)),double(x(itr)));
    end
    fprintf(fileID,'%s     %s %s\n','Execution time',num2str(execution),'seconds');
    fclose(fileID);
    
    
    
end

%%-----------------------------LU decomposition-------------------
if(find(selected_index==2))
        now1=tic();
    
        tol=str2double(get(handles.epsilon,'String'));
        maxIteration=str2double(get(handles.max_iteration,'String'));
        
        if(isempty(get(handles.max_iteration,'string')))
            maxIteration = 50;
        elseif(isempty(get(handles.epsilon,'string')))
            tol = 0.00001;
        end
   
NA = size(arr_final,1);
%concatinate A with identity matrix
AP = [arr_final eye(NA)]; 
for k = 1:NA - 1
    %Partial Pivoting at AP(k,k)
    %akx get greatest value in the row
    %kx get the number of row
    [akx, kx] = max(abs(AP(k:NA,k)));
    if akx < tol
        errordlg('the pivots less than tolerence','error')
        return;
    end
    mx = k+kx-1;
    %if the pivot not in the fist row then swap the 2 rows
    if kx > 1 
    tmp_row = AP(k,:);
    AP(k,:) = AP(mx,:);
    AP(mx,:) = tmp_row;
    end
    % LU decomposition
    %save both lower and upper in the same array
        for m = k + 1: NA
            if(AP(k,k) == 0)
            errordlg('The pivot cant be zero','error');
            return;
            end

        AP(m,k) = AP(m,k)/AP(k,k); 
        AP(m,k+1:NA) = AP(m,k + 1:NA)-AP(m,k)*AP(k,k + 1:NA); 
        end
end

%extract L and U matrixes from AP array
for m = 1:NA
    for n = 1:NA
            if m == n
                L(m,m) = 1;
                U(m,m) = AP(m,m);
            elseif m > n
                L(m,n) = AP(m,n);
                U(m,n) = 0;
            else
                L(m,n) = 0;
                U(m,n) = AP(m,n);
            end
    end
end

L
U

%FORWARD SUBSTITUTION
%forward substitution for a lower-triangular matrix equation Lx = B
N = size(L,1);  % ans = 3
if(L(1,1) == 0)
    errordlg('Singular matrix and No LU decomposition','error')
    return;
end
x1(1,:) = result_arr(1,:)/L(1,1);
for m = 2:N
x1(m,:) = (result_arr(m,:)-L(m,1:m - 1)*x1(1:m-1,:))/L(m,m);
end


%BACKWARD SUBSTITUTION
%backward substitution for a upper-triangular matrix equation Ux = B
N = size(U,2); % ans = 3
if(U(N,N) == 0)
    errordlg('Singular matrix and No LU decomposition','error')
    return;
end
x(N,:) = x1(N,:)/U(N,N);
for m = N-1: -1:1
x(m,:) = (x1(m,:) - U(m,m + 1:N)*x(m + 1:N,:))/U(m,m);
end

for i=1: len
        if isnan(x(i))
            errordlg('Singular matrix and No LU decomposition','Error');
            return;     
        elseif isinf(double(x(i)))
            errordlg('Singular matrix and No LU decomposition','Error');
            return;
        end
        
end

    set(handles.answer_box2,'String',double(x))
    execution=toc(now1);
    set(handles.time2,'string',[num2str(execution) ' seconds'])
    
    fileID = fopen('output.txt','a');
    fprintf(fileID,'%s\n','LU decomposition roots');
    fprintf(fileID,'%s     %s\n','var','value');
    for itr=1:len
        fprintf(fileID,'%s     %f\n',char(char_arr(itr)),double(x(itr)));
    end
    fprintf(fileID,'%s     %s %s\n','Execution time',num2str(execution),'seconds');
    fclose(fileID);
       
end

%%-----------------------------Gaussian-Jordan-------------------
if(find(selected_index==3))
    now1=tic();
    n= arr_size(1);
    matrix=arr_final(1:n,1:n);
    
    B = result_arr(1:n,:);
    B
    [m n]=size(matrix)
    
    nb=n+1
    matrix = [matrix B]
    %convert element below the major diagonal to zero
    for k=1:n-1
        for i=k+1:n
            m=matrix(i,k)/matrix(k,k)  %compute the kth column of m
            matrix(i,k:nb)=matrix(i,k:nb)-m*matrix(k,k:nb)    %compute An=Mn*An-1
        end
    end
    %convert element upper the major diagonal to zero
    for k=n:-1:2
        for i=k-1:-1:1
            m=matrix(i,k)/matrix(k,k)  %compute the kth column of m
            matrix(i,k:nb)=matrix(i,k:nb)-m*matrix(k,k:nb)    %compute An=Mn*An-1
        end
    end
    x=zeros(n,1);
    x(n)=matrix(n,nb)/matrix(n,n);
    %backward elimination
    for i=1:n
        matrix(i,:)= matrix(i,:)/matrix(i,i);
        
        if isnan(x(i))
            errordlg('INFINITE NO. OF SOLUTIONS','Error');
            
            return;
        elseif isinf(x(i))
            errordlg('SOLUTION DOES NOT EXIST!!','Error');
            return;
        else x(i)=matrix(i,n+1);
        end
    end
    set(handles.answer_box3,'String',double(x))
    execution=toc(now1);
    set(handles.time3,'string',[num2str(execution) ' seconds'])
    matrix
    
    fileID = fopen('output.txt','a');
    fprintf(fileID,'%s\n','Gaussian-Jordan roots');
    fprintf(fileID,'%s     %s\n','var','value');
    for itr=1:len
        fprintf(fileID,'%s     %f\n',char(char_arr(itr)),double(x(itr)));
    end
    fprintf(fileID,'%s     %s %s\n','Execution time',num2str(execution),'seconds');
    fclose(fileID);
end 
    %%-----------------------------Gauss-Seidel-------------------
    if(find(selected_index==4))
        
        now1=tic();
        A=arr_final;
        b=result_arr
        x=str2double(get(handles.initial_box,'String'));
        x_len=length(x)
        
        
        tol=str2double(get(handles.epsilon,'String'));
        maxIteration=str2double(get(handles.max_iteration,'String'));
        
        if(isempty(get(handles.max_iteration,'string')))
            maxIteration = 50;
        elseif(isempty(get(handles.epsilon,'string')))
            tol = 0.00001;
        end
        
        if(isempty(get(handles.initial_box,'string')))
            x=[zeros(len,1)]'
        elseif x_len ~= len
            errordlg('Missing input','Error');
        end
        
        isdom = true;
        dom = false;
        for r = 1:len
            rowdom = 2 * abs(A(r,r)) >= sum(abs(A(r,:)))
            isdom = isdom && rowdom;
            if 2 * abs(A(r,r)) > sum(abs(A(r,:)))
                diag(r)=1;
            end
        end
        
        if find(diag,1)
            dom = true;
        end
        isdom = isdom && dom;
        
        if ~isdom
            errordlg('Cannot use gauss seidel, the matrix is not strictly diagonally dominant','Error');
            return
        end
        
        itr=0;
        normVal=Inf;
        mm=1;
        nn=1;
        while normVal>tol && itr<maxIteration
            x_old=x;
            
            for i=1:len
                
                sigma=0;
                
                for j=1:i-1
                    sigma=sigma+A(i,j)*x(j);
                end
                
                for j=i+1:len
                    sigma=sigma+A(i,j)*x_old(j);
                end
                
                if A(i,i)==0
                    errordlg('Error cannot dvide by zero','Error');
                    return
                end
                
                x(i)=(1/A(i,i))*(b(i)-sigma)
                g(nn,mm)= x(i)
                mm=mm+1;
            end
            nn=nn+1;
            mm=1;
            itr=itr+1
            normVal=norm(x_old-x);
        end
        [rr cc]=size(g)
        %             itr=itr-1;
        for ii=1:cc
            axes(handles.graph);
            %         plot(itr,x(i));
            plot(g(:,ii))
            grid on;hold on;
            xlabel('iterations')
            ylabel('Roots')
        end
        
        set(handles.itable,'Data',double(g));
        set(handles.answer_box4,'String',double(x))
        set(handles.iterations,'String',itr)
        set(handles.precision,'String',normVal)
        execution=toc(now1);
        set(handles.time4,'string',[num2str(execution) ' seconds'])

        fileID = fopen('output.txt','a');
        fprintf(fileID,'%s\n','Gauss-Seidel roots');
        fprintf(fileID,'%s     %s\n','var','value');
        for i=1:len
            fprintf(fileID,'%s     %f\n',char(t(i)),double(x(i)));
        end
        fprintf(fileID,'%s    %f\n','Number of iterations',itr);
        fprintf(fileID,'%s    %f\n','Precision',normVal);
        
        for i=1:len
        fprintf(fileID,'%s         ',char(char_arr(i)));
        end
         fprintf(fileID,'\n');
        for ii=1:rr;
            for jj=1:cc
            fprintf(fileID,'%f  ',double(g(ii,jj)));
            
            end
            fprintf(fileID,'\n');
        end
        fprintf(fileID,'%s     %s %s\n','Execution time',num2str(execution),'seconds');
        fclose(fileID);


    end


% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in methods.
function methods_Callback(hObject, eventdata, handles)
% hObject    handle to methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methods contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methods


% --- Executes during object creation, after setting all properties.
function methods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in eq_table.
function eq_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to eq_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in eq_table.
function eq_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to eq_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in eq_box.
function eq_box_Callback(hObject, eventdata, handles)
% hObject    handle to eq_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns eq_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eq_box


% --- Executes during object creation, after setting all properties.
function eq_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eq_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in initial_box.
function initial_box_Callback(hObject, eventdata, handles)
% hObject    handle to initial_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns initial_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from initial_box


% --- Executes during object creation, after setting all properties.
function initial_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function result_Callback(hObject, eventdata, handles)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of result as text
%        str2double(get(hObject,'String')) returns contents of result as a double


% --- Executes during object creation, after setting all properties.
function result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function init_Callback(hObject, eventdata, handles)
% hObject    handle to init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of init as text
%        str2double(get(hObject,'String')) returns contents of init as a double


% --- Executes during object creation, after setting all properties.
function init_CreateFcn(hObject, eventdata, handles)
% hObject    handle to init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in add_init.
function add_init_Callback(hObject, eventdata, handles)
init = get(handles.init,'String');
row ={init};
oldData = get(handles.initial_box,'String');
newData=[oldData; row];
set(handles.initial_box,'String',newData);
set(handles.init,'String','');

% hObject    handle to add_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in answer_box2.
function answer_box2_Callback(hObject, eventdata, handles)
% hObject    handle to answer_box2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns answer_box2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from answer_box2


% --- Executes during object creation, after setting all properties.
function answer_box2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to answer_box2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in answer_box3.
function answer_box3_Callback(hObject, eventdata, handles)
% hObject    handle to answer_box3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns answer_box3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from answer_box3


% --- Executes during object creation, after setting all properties.
function answer_box3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to answer_box3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in answer_box1.
function answer_box1_Callback(hObject, eventdata, handles)
% hObject    handle to answer_box1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns answer_box1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from answer_box1


% --- Executes during object creation, after setting all properties.
function answer_box1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to answer_box1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in answer_box4.
function answer_box4_Callback(hObject, eventdata, handles)
% hObject    handle to answer_box4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns answer_box4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from answer_box4


% --- Executes during object creation, after setting all properties.
function answer_box4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to answer_box4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
