function Cr = pairCorr(image)

domains = imread('im123.png'); % read image

[rows, columns, numberOfColorChannels] = size(domains);
if numberOfColorChannels > 1
  % It's not really gray scale like we expected - it's color.
  % Convert it to gray scale by taking only the green channel.
  domains = domains(:, :, 3);   % Take green channel.
end

domains(domains>0) = 1;       % make sure its binary by setting 1 to values > 0
pp=sum(sum(domains));            % sum of total occupied sites
n = length(domains(:, 1));    % image size
rho=pp/(n*n)                 % density of system

for i=1:n
    line = domains(:, i);     % take one line...
    for j=1:n                 % and for each distance...
        s = line(1:end-n+j);
        Cr(i, j) = mean(s);   %...calculate Cr as mean  
    end
end

Cr = mean(Cr)/rho;

 figure(1); plot(Cr, '-o','LineWidth', 2);
 saveas(gcf,'plot-radf1.png');