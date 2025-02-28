clear; close all; tic

setenv('PATH', [getenv('PATH') ':/usr/local/Cellar/imagemagick/7.1.0-2_1/bin']);

%for linux "/usr/bin"
%for mac "/usr/local/Cellar/imagemagick/7.1.0-2_1/bin"
%myDir = uigetdir; %gets directory
%cd ('/Users/rakeshchatterjee/codes_/');
%cd ('/home/np1/staff/chatterjeerak/codes_desk/');

%==================
ens   = 5;
start = 3;
%==================


%---------------------- to define filenumber ---------------------------
p=1;
%===========
%filename1 = ['sh dospcon1-', num2str(p,'%04d'),'.sh']; system(filename1);
filename2 = ['sh dospcon2-', num2str(p,'%04d'),'.sh']; system(filename2);
%===========
myFiles1 = dir(fullfile('matimage*.png'));
M = length(myFiles1);
system('rm matimage*.*');
clear myFiles1;
cc_f = zeros(1,M+1-start);
a    = [1:M+1-start];
%-----------------------------------------------------------------------

%----------------- Average loop starts ---------------------------------
for p = 1:1:ens
%=========
%filename1 = ['sh dospcon1-', num2str(p,'%04d'),'.sh']; system(filename1);
filename2 = ['sh dospcon2-', num2str(p,'%04d'),'.sh']; system(filename2);
%=========

myFiles1 = dir(fullfile('matimage*.png'));
filenumber=0;


frame=0;
for t= 1:1:M
    frame = frame + 1;
    
  baseFileName = myFiles1(t).name;
  %fullFileName = fullfile(myDir, baseFileName);
  fullFileName = fullfile(baseFileName);
  %fprintf(1, 'Reading %s\n', baseFileName);
  num = importdata(fullFileName);   %or readtable
  filenumber = filenumber+1;
 

folder = fileparts(which(baseFileName));
fullFileName = fullfile(folder, baseFileName);

grayImage = imread(fullFileName);

[rows, columns, numberOfColorChannels] = size(grayImage);
if numberOfColorChannels > 1
  grayImage = grayImage(:, :, 3); % Take green channel.
end

grayImage(grayImage>0)=1;
grayImage = double(grayImage);
gimage0(:,:,frame) = grayImage;

end

gimage1 = gimage0(:,:,start:M);

Ai=zeros(columns);
for k=1:M+1-start
    Ai(:,:) = Ai(:,:) + gimage1(:,:,k);
end

Ai = Ai./(M+1-start);


fprintf('ens-> %d\n',p);


%------------------------------
%Calculation of Pearson correlation coefficient:
for t1=1:1:M+1-start
    cc1(t1) = corr2( ( gimage1(:,:,1) -Ai),( gimage1(:,:,t1) -Ai) );
    %cc1(t1) = corr2( ( (gimage1(:,:,1)/std2(gimage1(:,:,1))) -Ai/std2(Ai)),( (gimage1(:,:,t1)/std2(gimage1(:,:,t1))) -Ai/std2(Ai)) );
end

cc_f = cc_f + cc1;

system('rm matimage*.*');
clear myFiles1 baseFileName grayImage gimage0 gimage1 cc1;

end

cc_f = cc_f/ens;

%============
%system('rm dospcon1*.sh');
system('rm dospcon2*.sh');
%============

plot(cc_f,'-','LineWidth',2);
%xlim([1 lim]);
xlabel('t','FontSize',25);
ylabel('Pearson coeff','FontSize',25);

 data1= horzcat(a',cc_f');
%============= 
        saveas(gcf,'plotL102bnew.png');
 writematrix(data1,'dataL102bnew.txt','Delimiter','tab');
%=============
toc
