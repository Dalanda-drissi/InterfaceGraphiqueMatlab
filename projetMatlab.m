function varargout = projetMatlab(varargin)
% PROJETMATLAB MATLAB code for projetMatlab.fig
%      PROJETMATLAB, by itself, creates a new PROJETMATLAB or raises the existing
%      singleton*.
%
%      H = PROJETMATLAB returns the handle to a new PROJETMATLAB or the handle to
%      the existing singleton*.
%
%      PROJETMATLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJETMATLAB.M with the given input arguments.
%
%      PROJETMATLAB('Property','Value',...) creates a new PROJETMATLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before projetMatlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to projetMatlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help projetMatlab

% Last Modified by GUIDE v2.5 24-Dec-2020 14:11:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @projetMatlab_OpeningFcn, ...
                   'gui_OutputFcn',  @projetMatlab_OutputFcn, ...
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


% --- Executes just before projetMatlab is made visible.
function projetMatlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to projetMatlab (see VARARGIN)
set(handles.EtOuNotPanel,'visible','off');
set(handles.axes10,'visible','off');
  set(handles.axes11,'visible','off');
  set(handles.axes3,'visible','off');
  set(handles.Average,'visible','off');
  set(handles.Gaussian,'visible','off');
  set(handles.pnsr,'visible','off');
   set(handles.bruitTxt,'visible','off');
   set(handles.psnrTXT,'Visible','off');
   set(handles.psnrTXT1,'Visible','off');

  set(handles.CodeSourcePanel,'visible','off');
 set(handles.ImageOriginal,'visible','off');
  set(handles.bruitTxt,'visible','off');
  set(handles.ImageModif,'visible','off');
  set(handles.HistoTxt,'visible','off');
% Choose default command line output for projetMatlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projetMatlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = projetMatlab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
function ResetAxes(handles)
     cla(handles.axes3,'reset');    
     cla(handles.axes10,'reset');
     cla(handles.axes11,'reset');
     cla(handles.axes13,'reset');
     cla(handles.axes14,'reset');
     cla(handles.axes15,'reset');
     cla(handles.axes16,'reset');
     cla(handles.axes18,'reset');
     cla(handles.axes19,'reset');
function setOFF (handles)
 set(handles.CodeSourcePanel,'visible','off');
  set(handles.axes10,'visible','off');
  set(handles.axes11,'visible','off');
  set(handles.Average,'visible','off');
  set(handles.Gaussian,'visible','off');
  set(handles.pnsr,'visible','off');
   set(handles.bruitTxt,'visible','off');
   set(handles.psnrTXT,'Visible','off');
set(handles.psnrTXT1,'Visible','off');
set(handles.ImageOriginal,'visible','off');
  set(handles.bruitTxt,'visible','off');
  set(handles.ImageModif,'visible','off');
  set(handles.HistoTxt,'visible','off');
function setONN (handles)
  set(handles.axes10,'visible','on');
  set(handles.axes11,'visible','on');
 
  set(handles.Average,'visible','on');
  set(handles.Gaussian,'visible','on');
  set(handles.pnsr,'visible','on');
  set(handles.bruitTxt,'visible','on');
  
% --- Executes on button press in chooseImagek.
function chooseImagek_Callback(hObject, eventdata, handles)
ResetAxes(handles);
set(handles.axes3,'visible','off');
setOFF (handles);
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
global selectedImage
[filename,filepath]=uigetfile({'*.bmp'},'Select and image');
selectedImage = imread(strcat(filepath, filename));
axes(handles.axes1);
imshow(selectedImage);


% --- Executes on button press in RGBtoGRAY.
function RGBtoGRAY_Callback(hObject, eventdata, handles)
ResetAxes(handles);
setOFF (handles);

set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
 global selectedImage
   set(handles.Panel,'Title','RGB to GRAY');
    GrayImage=rgb2gray(selectedImage);
    axes(handles.axes3);
    imshow(GrayImage);
 
    


% --- RGB to Black/white.
function RGBTOBW_Callback(hObject, eventdata, handles)
ResetAxes(handles);
    setOFF (handles);
set(handles.axes10,'visible','off');
  set(handles.axes11,'visible','off');
  
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
global selectedImage
set(handles.Panel,'Title','RGB to Black / White');
WhiteBlackImg = im2bw(selectedImage);
axes(handles.axes3);
imshow(WhiteBlackImg);

% --- Executes on button press in rgb2hsv.
function rgb2hsv_Callback(hObject, eventdata, handles)

ResetAxes(handles);
setOFF (handles);
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
global selectedImage
set(handles.Panel,'Title','RGB to HSV');
ImageHsv= rgb2hsv(selectedImage);
axes(handles.axes3);
imshow (ImageHsv);

% --- Executes on button press in Assombrir.
function Assombrir_Callback(hObject, eventdata, handles)
ResetAxes(handles);
    setOFF (handles);
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
global selectedImage
set(handles.Panel,'Title','Assombrir une image');
c= sqrt(double(selectedImage));
axes(handles.axes3);
imshow(uint8 (c))

% --- Executes on button press in Eclairir.
function Eclairir_Callback(hObject, eventdata, handles)
setOFF (handles);
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
global selectedImage
set(handles.Panel,'Title','Eclairir une image');
axes(handles.axes3);
imshow(uint8 (selectedImage.^2))

% --- Executes on button press in Filtrage.
function Filtrage_Callback(hObject, eventdata, handles)
ResetAxes(handles);
set(handles.CodeSourcePanel,'visible','off');
set(handles.EtOuNotPanel,'visible','off');
set(handles.Panel,'visible','on');
global selectedImage
global B
setONN (handles);
set(handles.Panel,'Title','Filtrage');
[filename,filepath]=uigetfile({'*.bmp'},'Select and image');
selectedImage = imread(strcat(filepath, filename));
axes(handles.axes1);
imshow(selectedImage);
B=imnoise(selectedImage,'salt & pepper');%ajouter  bruit
axes(handles.axes3);
imshow(B);

% --- Executes on button press in CodeSource.
function CodeSource_Callback(hObject, eventdata, handles)

    str = get (handles.Panel,'Title');
  
    
    if (strncmp(str,'RGB to GRAY',10))
        text=sprintf('axes(handles.axes3); \n imshow(uint8 (selectedImage^.2));');
        set(handles.CodeSourcePanel,'visible','on');
        set(handles.text,'String',text);
    end
    if (strncmp(str,'RGB to Black / White',19))
         text=sprintf('WhiteBlackImg = im2bw(selectedImage); \n axes(handles.axes3); \n imshow(WhiteBlackImg);');
        set(handles.CodeSourcePanel,'visible','on');
        set(handles.text,'String',text);
    end
    if (strncmp(str,'RGB to HSV',10))
         text=sprintf('ImageHsv= rgb2hsv(selectedImage);\n axes(handles.axes3);\n imshow (ImageHsv);');
        set(handles.CodeSourcePanel,'visible','on');
        set(handles.text,'String',text);
    end
    if (strncmp(str,'Assombrir une image',19))
         text=sprintf('c= sqrt(double(selectedImage)); \n axes(handles.axes3);\n imshow(uint8 (c))');
        set(handles.CodeSourcePanel,'visible','on');
        set(handles.text,'String',text);
    end
    if (strncmp(str,'Eclairir une image',19))
         text=sprintf('c= sqrt(double(selectedImage)); \n axes(handles.axes3);\n imshow(uint8 (c))');
        set(handles.CodeSourcePanel,'visible','on');
        set(handles.text,'String',text);
    end
     if (strncmp(str,'Filtrage/Average',15))
         text=sprintf('B=imnoise(selectedImage,''salt & pepper''); \n h1=fspecial(''average'',[3 3]);\n h2=fspecial(''average'',[5 5]);\n h3=fspecial(''average'',[7 7]);\n A1=imfilter(B,h1);\n  A2=imfilter(B,h2);\n A3=imfilter(B,h3);\n subplot(231);\n imshow(A);\n title(''Image'');\n subplot(232);\n imshow(B);\n title(''ajouter  bruit'');\n subplot(234);\n imshow(A1);\n subplot(235);\n imshow(A2);\n subplot(236);\n imshow(A3);');
        set(handles.CodeSourcePanel,'visible','on');
        set(handles.text,'String',text);
     end
    if (strncmp(str,'Filtrage/Gaussian',15))
         text=sprintf('B=imnoise(selectedImage,''salt & pepper'');\n h1=fspecial(''gaussian'',[5 5], 0.9); (gamma ,  taux  qui controle le bruit entre [0, 1]) \n A1=imfilter(B,h1);\n axes(handles.axes10);\n imshow(A1);');
        set(handles.CodeSourcePanel,'visible','on');
        set(handles.text,'String',text);
    end
     if (strncmp(str,'Filtrage/Pnsr',15))
         text=sprintf('function psnr = psnr( X,Y )\n [a,b]=size(X);\n somme =0;\n z=double(X)-double(Y);\n for i=1:a, \n  for j=1:b,  \n    somme=somme+z(i,j)*z(i,j); \n end\n end\n er=somme/(a*b);\n psnr=(10* log10((255*255/er)));'); 
         set(handles.CodeSourcePanel,'visible','on'); 
         set(handles.text,'String',text);
    end
     if (strncmp(str,'Filtrage/Image Colorée',15))
         text=sprintf('A=imread(''peppers.png'');\n B=imnoise(A,''salt & pepper'');\n r=B(:,:,1);\nv=B(:,:,2);\nb= B(:,:,3);\nh=fspecial(''average'',[3 3]);\nr1=imfilter(r,h);\n v1=imfilter(v,h);\nb1=imfilter(b,h);\nkfilt(:,:,1) = r1;\nkfilt(:,:,2) = v1;\nkfilt(:,:,3) = b1;\nsubplot(131);imshow(A)\nsubplot(132);imshow(B)\nsubplot(133);imshow(kfilt);'); 
         set(handles.CodeSourcePanel,'visible','on'); 
         set(handles.text,'String',text);
    end
    
     if (strncmp(str,'Rehaussement d’images/Inversion',15))
         text=sprintf('A = imread(''image/scan.bmp'');\nsubplot(231); imshow(A);title(''image originale'');\nsubplot(232) ; imhist(A) ; title(''son histogramme'')\nImage=imread(''image/scan.bmp'');\n[L,l]=size(Image);\nfor i=1:L\n for j=1:l\n J1(i,j)=255-Image(i,j);\n end\nend\nsubplot(234) ; imshow(J1) ; title(''image inversÃ©e'')\nsubplot(235) ; imhist(J1) ; title(''son histogramme'')'); 
         set(handles.CodeSourcePanel,'visible','on'); 
         set(handles.text,'String',text);
     end
    if (strncmp(str,'Rehaussement d’images/Egalisation',25))
         text=sprintf('\nA = imread(''tire.tif'');\nsubplot(231); imshow(A);title(''image originale'');\nsubplot(232) ; imhist(A) ; title(''son histogramme'')\nA1 = histeq(A); \nsubplot(234) ; imshow(A1) ; title('' image egalisÃ©e'')\nsubplot(235) ; imhist(A1) ; title(''son histogramme'')'); 
         set(handles.CodeSourcePanel,'visible','on'); 
         set(handles.text,'String',text);
    end
    if (strncmp(str,'Rehaussement d’images/Recadrage',25))
         text=sprintf('\nA = imread(''pout.tif'');\nsubplot(231); imshow(A);title(''image originale'');\nsubplot(232) ; imhist(A) ; title(''son histogramme'')\nA1 = imadjust(A);\nsubplot(234) ; imshow(A1) ; title(''image recadrÃ©e'')\nsubplot(235) ; imhist(A1) ; title(''son histogramme'')'); 
         set(handles.CodeSourcePanel,'visible','on'); 
         set(handles.text,'String',text);
    end
    




 
function psnr = psnr( X,Y )

[a,b]=size(X);
somme =0;
z=double(X)-double(Y);
for i=1:a,
   for j=1:b, 
       somme=somme+z(i,j)*z(i,j);
   end
end
er=somme/(a*b);
psnr=(10* log10((255*255/er)));



% --- Executes on button press in Average.
function Average_Callback(hObject, eventdata, handles)
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
set(handles.Panel,'Title','Filtrage/Average');
global B
h1=fspecial('average',[3 3]);
A1=imfilter(B,h1);
axes(handles.axes11);
imshow(A1);



% --- Executes on button press in Gaussian.
function Gaussian_Callback(hObject, eventdata, handles)
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
global B
set(handles.Panel,'Title','Filtrage/Gaussian');
h1=fspecial('gaussian',[5 5], 0.5); %gamma ,  taux  qui controle le bruit entre [0, 1]
A1=imfilter(B,h1);
axes(handles.axes10);
imshow(A1);


% --- Executes on button press in pnsr.
function pnsr_Callback(hObject, eventdata, handles)
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
set(handles.Panel,'Title','Filtrage/Pnsr');
set(handles.psnrTXT,'Visible','on');
set(handles.psnrTXT1,'Visible','on');
global selectedImage
global B
h1=fspecial('gaussian',[5 5], 0.5); %gamma ,  taux  qui controle le bruit entre [0, 1]
h2=fspecial('average',[3 3]);
A1=imfilter(B,h1);
A2=imfilter(B,h2);
psnr1=psnr(selectedImage,A1);
psnr2=psnr(selectedImage,A2);
text=sprintf(num2str(psnr1)); % Gaussian 
text1=sprintf(num2str(psnr2)); % Average
set(handles.psnrTXT,'String',text);
set(handles.psnrTXT1,'String',text1);


% --- Executes on button press in FiltrageImageColore.
function FiltrageImageColore_Callback(hObject, eventdata, handles)
    global selectedImage;
    ResetAxes(handles);
set(handles.Panel,'visible','on');
set(handles.EtOuNotPanel,'visible','off');
set(handles.Panel,'Title','Filtrage/Image Colorée');
setOFF (handles);
set(handles.axes11,'visible','on');
 cla(handles.axes10,'reset');
set(handles.axes10,'visible','off');

[filename,filepath]=uigetfile({'*.bmp'},'Select and image');
selectedImage = imread(strcat(filepath, filename));
axes(handles.axes1);
imshow(selectedImage);
B=imnoise(selectedImage,'salt & pepper');%ajouter  bruit
axes(handles.axes3);
imshow(B);
 r=B(:,:,1);
 v=B(:,:,2);
 b= B(:,:,3);
 h=fspecial('average',[3 3]);
 r1=imfilter(r,h);
 v1=imfilter(v,h);
 b1=imfilter(b,h);
 kfilt(:,:,1) = r1;
 kfilt(:,:,2) = v1;
 kfilt(:,:,3) = b1;
 axes(handles.axes11);
imshow(kfilt);


% --- Executes on button press in EtOuNot.
function EtOuNot_Callback(hObject, eventdata, handles)
 ResetAxes(handles);
set(handles.EtOuNotPanel,'visible','on');
set(handles.CodeSourcePanel2,'visible','off'); 
set(handles.Panel,'visible','off');


% --- Executes on button press in NotBtn.
function NotBtn_Callback(hObject, eventdata, handles)
    
global H;
set(handles.EtOuNotPanel,'Title','Not');
set(handles.Panel,'visible','off');
[m,n]= size(H);
for i = 1:m
    for j=1:n
        if (H(i,j)==1 )
           
            F(i,j)=0;
        else
            F(i,j)=1;
        end
    end
end
 axes(handles.axes19);
imshow(F);


% --- Executes on button press in OuBtn.
function OuBtn_Callback(hObject, eventdata, handles)
   
global H;
global K;
global L;
set(handles.EtOuNotPanel,'Title','Ou');
set(handles.Panel,'visible','off');
[m,n]= size(H);
for i = 1:m
    for j=1:n
        if (H(i,j)==1 || K(i,j) == 1 || L(i,j)==1)
            M(i,j)=1;
        else
            M(i,j)=0;
        end
    end
end
axes(handles.axes15);
imshow(M);
% --- Executes on button press in EtBtn.
function EtBtn_Callback(hObject, eventdata, handles)
  
global H;
global K;
global L;
set(handles.EtOuNotPanel,'Title','Et');
set(handles.Panel,'visible','off');
[m,n]= size(H);
for i = 1:m
    for j=1:n
        if (H(i,j)==1 && K(i,j) == 1 && L(i,j)==1)
            G(i,j)=1;
        else
            G(i,j)=0;
        end
    end
end
axes(handles.axes16);
imshow(G);
% --- Executes on button press in Image1.
function Image1_Callback(hObject, eventdata, handles)
global H;
set(handles.Panel,'visible','off');
[filename,filepath]=uigetfile({'*.bmp'},'Select and image');
H = imread(strcat(filepath, filename));
axes(handles.axes13);
imshow(H);


% --- Executes on button press in Image2.
function Image2_Callback(hObject, eventdata, handles)
global K;
set(handles.Panel,'visible','off');
[filename,filepath]=uigetfile({'*.bmp'},'Select and image');
K = imread(strcat(filepath, filename));
axes(handles.axes14);
imshow(K);

% --- Executes on button press in Image3.
function Image3_Callback(hObject, eventdata, handles)
global L;
set(handles.Panel,'visible','off');
set(handles.Panel,'visible','off');
[filename,filepath]=uigetfile({'*.bmp'},'Select and image');
L = imread(strcat(filepath, filename));
axes(handles.axes18);
imshow(L);


% --- Executes on button press in CodeSourceEtOuNot.
function CodeSourceEtOuNot_Callback(hObject, eventdata, handles)
str2 = get (handles.EtOuNotPanel,'Title');
    if (strncmp(str2,'Et',2))
         text=sprintf('[m,n]= size(H);\nfor i = 1:m \nfor j=1:n  \n  if (H(i,j)==1 && K(i,j) == 1 && L(i,j)==1)    \n  G(i,j)=1; \n else \n   G(i,j)=0; \n   end \n  end\nend'); 
         set(handles.CodeSourcePanel2,'visible','on'); 
         set(handles.TextCodeSourceEON,'String',text);
    end
    if (strncmp(str2,'Ou',2))
         text=sprintf('for i = 1:m\n for j=1:n\n  if (H(i,j)==1 || K(i,j) == 1 || L(i,j)==1) \n  M(i,j)=1;\n  else\nM(i,j)=0;\n  end\n end\nend'); 
         set(handles.CodeSourcePanel2,'visible','on'); 
         set(handles.TextCodeSourceEON,'String',text);
    end
    if (strncmp(str2,'Not',2))
         text=sprintf('[m,n]= size(H);\n for i = 1:m\n  for j=1:n\nif (H(i,j)==1 )\n  F(i,j)=0;\n else\n   F(i,j)=1;\n end\nend\nend'); 
         set(handles.CodeSourcePanel2,'visible','on'); 
         set(handles.TextCodeSourceEON,'String',text);
    end


% --- Executes on button press in Inversion.
function Inversion_Callback(hObject, eventdata, handles)
     ResetAxes(handles);
global selectedImage
set(handles.Panel,'Title','Rehaussement d’images/Inversion');
 RehaussementImageHandles(handles);
  set(handles.ImageModif,'String','Image Inverse');
  axes(handles.axes1);
  imshow(selectedImage);
  axes(handles.axes3);
  imhist(selectedImage);
  [L,l]=size(selectedImage);
for i=1:L
   for j=1:l      
      J1(i,j)=255-selectedImage(i,j);            
   end
end
axes(handles.axes11);
imshow(J1);
axes(handles.axes10);
imhist(J1);


% --- Executes on button press in Egalisation.
function Egalisation_Callback(hObject, eventdata, handles)
     ResetAxes(handles);
global selectedImage
set(handles.Panel,'Title','Rehaussement d’images/Egalisation');
  RehaussementImageHandles(handles);
  set(handles.ImageModif,'String','Image Egalisee');
  axes(handles.axes1);
  imshow(selectedImage);
  axes(handles.axes3);
  imhist(selectedImage);
  A1 = histeq(selectedImage); 
  axes(handles.axes11);
  imshow(A1);
  axes(handles.axes10);
  imhist(A1);
  
% --- Executes on button press in Recadrage.
function Recadrage_Callback(hObject, eventdata, handles)
     ResetAxes(handles);
global selectedImage
set(handles.Panel,'Title','Rehaussement d’images/Recadrage');
  RehaussementImageHandles(handles);
  set(handles.ImageModif,'String','Image Recadrée');
  axes(handles.axes1);
  imshow(selectedImage);
  axes(handles.axes3);
  imhist(selectedImage);
  A1 = imadjust(selectedImage); 
  axes(handles.axes11);
  imshow(A1);
  axes(handles.axes10);
  imhist(A1);
  
  function RehaussementImageHandles(handles)
      set(handles.CodeSourcePanel,'visible','off'); 
  set(handles.axes10,'visible','on');
  set(handles.axes11,'visible','on');
  set(handles.Average,'visible','off');
  set(handles.Gaussian,'visible','off');
  set(handles.pnsr,'visible','off');
  set(handles.ImageOriginal,'visible','on');
  set(handles.bruitTxt,'visible','on');
  set(handles.ImageModif,'visible','on');
  set(handles.HistoTxt,'visible','on');
  set(handles.EtOuNotPanel,'visible','off');
  set(handles.Panel,'visible','on');
  set(handles.ImageOriginal,'String','Image Original');
  set(handles.bruitTxt,'String','son histogramme'); 
