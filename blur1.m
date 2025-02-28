clear; close all; tic
setenv('PATH', [getenv('PATH') ':/usr/local/Cellar/imagemagick/7.1.0-2_1/bin']);

%===================
A = 'mat1imagel102-002000.png';
B = 'mat2imagel102-002000.png';

sigma   = 5.0;
scale   = 4.5;
channel = 3;
%===================

[I21,cmap1] = imread(A);
I21= ind2rgb(I21,cmap1);
I3a=imgaussfilt(I21,sigma);
I3a=I3a*scale;

%figure(1); imshow(I3a); 
imwrite(I3a, "filter-g2.png", "png");

%I4 = imadjust(I3,[.1 .2 0; 0.2 0.3 1],[]);figure(2); imshow(I4);


[I22,cmap2] = imread(B);
I22= ind2rgb(I22,cmap2);
I3b=imgaussfilt(I22,sigma);
I3b=I3b*scale;

%figure(2); imshow(I3b); 
imwrite(I3b, "filter-r2.png", "png");

%I4 = imadjust(I3,[.1 .2 0; 0.2 0.3 1],[]);figure(2); imshow(I4);

%----------- Merging images ----------
%C2 = I3a + I3b; 
C2 = imfuse(I3a,I3b,'falsecolor','ColorChannels',[2 1 0]);
%C2 = imfuse(I3a,I3b,'ColorChannels',[2 1 0]);
%C2=C2*scale;
figure(3); imshow(C2);
imwrite(C2, "fusedimg.png", "png");
%--------------------------------------

%---- Measuring Correlation coefficient ----
A1=imread('filter-g2.png');
A1 = A1(:, :, channel);
B1=imread('filter-r2.png');
B1 = B1(:, :, channel);


coef = corr2(A1,B1)
fprintf('%c',A);fprintf('\n')

