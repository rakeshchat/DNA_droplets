clear; close all; tic
setenv('PATH', [getenv('PATH') ':/usr/local/Cellar/imagemagick/7.1.0-2_1/bin']);

%==================
sigma   = 5.0;
scale   = 8.0;
%==================
ens     = 1;
channel = 3;
M       = 1000;
tmax    = power(10,8);
data    = power(10,3);
gap     = tmax/data;
%==================


for p= 1:1:ens

    
%before this, execute the program split-convert.c
%filename = ['sh dospcon1-', num2str(p,'%04d'),'.sh'];


coef = zeros(1,M);
a=[1:gap:tmax];
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


filename_img1 = ['imageblur1_',num2str(t,'%06d'),'.png'];
imwrite(Gs1, filename_img1, "png");

filename_img2 = ['imageblur2_',num2str(t,'%06d'),'.png'];
imwrite(Gs2, filename_img2, "png");

A1=imread(filename_img1);  
B1=imread(filename_img2);

A1 = A1(:, :, channel);
B1 = B1(:, :, channel);

coef(t) = coef(t) + corr2(A1,B1);

end

system('rm mat1image*.png');
system('rm mat2image*.png');
system('rm imageblur*.png');


clear myFiles1 myFiles2 baseFileName1 baseFileName2 grayImage1 grayImage2;

end

data1= horzcat(a',(coef/ens)');   
writematrix(data1,'data_rrfixrg2.txt','Delimiter','tab');



toc
