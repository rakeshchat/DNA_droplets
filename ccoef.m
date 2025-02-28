clear; close all; tic

setenv('PATH', [getenv('PATH') ':/usr/local/Cellar/imagemagick/7.1.0-2_1/bin']);


%for linux "/usr/bin"
%for mac "/usr/local/Cellar/imagemagick/7.1.0-2_1/bin"

%myDir = uigetdir; %gets directory
%cd ('/Users/rakeshchatterjee/codes_/');
%cd ('/home/np1/staff/chatterjeerak/codes_desk/');


%system('./dosplit.sh');
%system('./convert.sh');

%==========================================
myFiles1 = dir(fullfile('matimage*.png'));
filenumber=0;
M = length(myFiles1);
a=[1:M];


frame=0;
for t= 1:1:M
    frame = frame + 1;
    
  baseFileName = myFiles1(t).name;
  %fullFileName = fullfile(myDir, baseFileName);
  fullFileName = fullfile(baseFileName);
  fprintf(1, 'Reading %s\n', baseFileName);
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

Ai=zeros(columns);
for k=1:M
    Ai(:,:) = Ai(:,:) + gimage0(:,:,k);
end

Ai = Ai./M;

%======================================

%Calculation of Pearson correlation coefficient:
for t1=1:1:M
    %cc1(t1) = corr2( ( (gimage0(:,:,2)/std2(gimage0(:,:,1))) -Ai),( (gimage0(:,:,t1)/std2(gimage0(:,:,t1))) -Ai) );
    
    cc1(t1) = corr2( ( (gimage0(:,:,1)/std2(gimage0(:,:,1))) -Ai),( (gimage0(:,:,t1)/std2(gimage0(:,:,t1))) -Ai) );
    
end


plot(cc1,'-','LineWidth',2);
%xlim([1 lim]);
xlabel('t','FontSize',25);
ylabel('Pearson coeff','FontSize',25);

 data1= horzcat(a',cc1');
 
        saveas(gcf,'plot1L100ab.png');
 writematrix(data1,'data1L100ab.txt','Delimiter','tab');

toc



%system('rm yimage*.png');
