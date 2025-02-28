clear; close all; tic

setenv('PATH', [getenv('PATH') ':/usr/local/Cellar/imagemagick/7.1.0-2_1/bin']);

%for linux "/usr/bin"
%for mac "/usr/local/Cellar/imagemagick/7.1.0-2_1/bin"
%myDir = uigetdir; %gets directory
%cd ('/Users/rakeshchatterjee/codes_/');
%cd ('/home/np1/staff/chatterjeerak/codes_desk/');

%==================
ens   = 2;
start = 1;
%==================


%---------------------- to define filenumber ---------------------------
p=1;
%===========
filename1 = ['sh dospcon1-', num2str(p,'%04d'),'.sh']; system(filename1);
%===========
myFiles1 = dir(fullfile('mat1image*.png'));
M = length(myFiles1);
system('rm mat1image*.*');
clear myFiles1;
cc_f = zeros(1,M+1-start);
a    = [1:M+1-start];
%-----------------------------------------------------------------------

%----------------- Average loop starts ---------------------------------
for p = 1:1:ens
%=========
filename1 = ['sh dospcon1-', num2str(p,'%04d'),'.sh']; system(filename1);
filename2 = ['sh dospcon2-', num2str(p,'%04d'),'.sh']; system(filename2);
%=========

myFiles1 = dir(fullfile('mat1image*.png'));
filenumber1=0;


frame=0;
for t= 1:1:M
    frame = frame + 1;
    
  baseFileName1 = myFiles1(t).name;
  %fullFileName1 = fullfile(myDir, baseFileName1);
  fullFileName1 = fullfile(baseFileName1);
  %fprintf(1, 'Reading %s\n', baseFileName1);
  num = importdata(fullFileName1);   %or readtable
  filenumber1 = filenumber1+1;
 

folder = fileparts(which(baseFileName1));
fullFileName1 = fullfile(folder, baseFileName1);

grayImage1 = imread(fullFileName1);

[rows, columns, numberOfColorChannels] = size(grayImage1);
if numberOfColorChannels > 1
  grayImage1 = grayImage1(:, :, 3); % Take green channel.
end

%imshow(double(grayImage1));

grayImage1(grayImage1>0)=1;
grayImage1 = double(grayImage1);
gimage01(:,:,frame) = grayImage1;

end

gimage1 = gimage01(:,:,start:M);

Ai=zeros(columns);
for k=1:M+1-start
    Ai(:,:) = Ai(:,:) + gimage1(:,:,k);
end

Ai = Ai./(M+1-start);

%---------------------
myFiles2 = dir(fullfile('mat2image*.png'));
filenumber2=0;


frame=0;
for t= 1:1:M
    frame = frame + 1;
    
  baseFileName2 = myFiles2(t).name;
  %fullFileName2 = fullfile(myDir, baseFileName2);
  fullFileName2 = fullfile(baseFileName2);
  %fprintf(1, 'Reading %s\n', baseFileName2);
  num = importdata(fullFileName2);   %or readtable
  filenumber2 = filenumber2+1;
 

folder = fileparts(which(baseFileName2));
fullFileName2 = fullfile(folder, baseFileName2);

grayImage2 = imread(fullFileName2);

[rows, columns, numberOfColorChannels] = size(grayImage1);
if numberOfColorChannels > 1
  grayImage2 = grayImage2(:, :, 3); % Take green channel.
end

%imshow(double(grayImage2));

grayImage2(grayImage2>0)=1;
grayImage2 = double(grayImage2);
gimage02(:,:,frame) = grayImage2;

end

gimage2 = gimage02(:,:,start:M);

Aj=zeros(columns);
for k=1:M+1-start
    Aj(:,:) = Aj(:,:) + gimage2(:,:,k);
end

Aj = Aj./(M+1-start);

%------------------------



fprintf('ens-> %d\n',p);


%------------------------------
%Calculation of Pearson correlation coefficient:
for t1=1:1:M+1-start
    
    cc1(t1) = corr2( ( gimage1(:,:,t1) -Ai),( gimage2(:,:,t1) -Aj) );

end

cc_f = cc_f + cc1;

system('rm mat*image*.*');
clear myFiles1 baseFileName1 grayImage1 gimage01 gimage1 myFiles2 baseFileName2 grayImage2 gimage02 gimage2 cc1;

end

cc_f = cc_f/ens;

%============
%system('rm dospcon1*.sh');
%system('rm dospcon2*.sh');
%============

plot(cc_f,'-','LineWidth',2);
%xlim([1 lim]);
xlabel('t','FontSize',25);
ylabel('Pearson coeff','FontSize',25);

 data1= horzcat(a',cc_f');
%============= 
        saveas(gcf,'plL52.png');
 writematrix(data1,'dtL52.txt','Delimiter','tab');
%=============
toc
