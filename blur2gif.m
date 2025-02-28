clear; close all; tic
setenv('PATH', [getenv('PATH') ':/usr/local/Cellar/imagemagick/7.1.0-2_1/bin']);

%==================
sigma = 5.0;
scale = 4.5;
M     = 1000;
%==================


myFiles1 = dir(fullfile('mat1image-*.png'));
myFiles2 = dir(fullfile('mat2image-*.png'));
filenumber1=0; filenumber2=0;


frame=0;
for t= 1:1:M
    frame = frame + 1;
    
  baseFileName1 = myFiles1(t).name;
  %fullFileName1 = fullfile(myDir, baseFileName1);
  fullFileName1 = fullfile(baseFileName1);
  %fprintf(1, 'Reading %s\n', baseFileName1);
  num1 = importdata(fullFileName1);   %or readtable
  filenumber1 = filenumber1+1;


  baseFileName2 = myFiles2(t).name;
  %fullFileName2 = fullfile(myDir, baseFileName2);
  fullFileName2 = fullfile(baseFileName2);
  %fprintf(1, 'Reading %s\n', baseFileName2);
  num2 = importdata(fullFileName2);   %or readtable
  filenumber2 = filenumber2+1;
 

folder = fileparts(which(baseFileName1));
fullFileName1 = fullfile(folder, baseFileName1);

folder = fileparts(which(baseFileName2));
fullFileName2 = fullfile(folder, baseFileName2);


[grayImage1,cmap] = imread(fullFileName1);
grayImage1= ind2rgb(grayImage1,cmap);
Gs1=imgaussfilt(grayImage1,sigma);
Gs1=Gs1*scale;

[grayImage2,cmap] = imread(fullFileName2);
grayImage2= ind2rgb(grayImage2,cmap);
Gs2=imgaussfilt(grayImage2,sigma);
Gs2=Gs2*scale;


%J1 = imadjust(Gs1,[.2 .3 0; .5 0.6 1],[]);
%J2 = imadjust(Gs2,[.2 .3 0; .9 1.0 1],[]);


C1 = imfuse(Gs1,Gs2,'falsecolor','ColorChannels',[2 1 0]);

filename_img = ['imageblur_',num2str(t,'%06d'),'.png'];

imwrite(C1, filename_img, "png");

end

system('convert -delay 14 -loop 1 imageblur_*.png gblurs5045.gif');
%system('rm imageblur_*.png');

toc
