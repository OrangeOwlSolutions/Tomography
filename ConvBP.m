%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVOLUTION BACKPROJECTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ConvBP(dataFileName, reconstructionFile, x, y)

fontSize = 18;

load(dataFileName);

% --- Smallest power of 2 greater than numDetectors
n1 = ceil(log(numDetectors) / log(2)); 
% N = 2^n1;                  
N = 2^(n1 + 1);

% --- Inter-detector distance
dt = (max(detectorPositions) - min(detectorPositions)) / (numDetectors - 1);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZERO-PADDED FFTs OF THE PROJECTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filteredProjectionSpectra = fft(sinogram.', N);        % --- Padded FFT is performed for every column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCTING THE BANDLIMITED FILTER %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Filter input response. N/2 instead of N is used for zero padding of h to avoid interperiod interference. 
h = zeros(N, 1);
h(1)                           = 1 / 4;                                            % --- n = 0
h(2 :  2 : (N / 2))            = -1 ./ (((2 : 2 : (N / 2)) - 1).^2 * pi^2);        % --- n odd
h(N : -2 : (N - N / 2 + 2))    = h(2 : 2 : (N / 2));                               

H = fft(h) / dt;                    % --- Padded FFT is performed for every column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTERING IN THE SPECTRAL DOMAIN %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filteredProjectionSpectra = H .* filteredProjectionSpectra;

filteredProjections = real(ifft(filteredProjectionSpectra))';   % --- IFFT is performed along every column

%%%%%%%%%%%%%%%%%%
% BACKPROJECTION %
%%%%%%%%%%%%%%%%%%
[XX, YY] = meshgrid(x, -y);
Reconstruction = zeros(size(XX));

for currentProjection = 1 : numProjections,     
    
    theta = angles(currentProjection);

    t1 = XX * cos(theta) + YY * sin(theta);

    currentReconstruction = filteredProjections(currentProjection, 1 : numDetectors);

    Reconstruction = Reconstruction + interp1(detectorPositions, currentReconstruction, t1, 'linear');   
end

%%%%%%%%%%%%%%%%%%%%%%
% GRAPHS AND SAVINGS %
%%%%%%%%%%%%%%%%%%%%%%
figure
imagesc(x, -y, Reconstruction / max(max(Reconstruction)), [0 1])
set(gca, 'YDir', 'Normal')
colormap(gray)
axis square
title('Convolution BackProjection', 'FontSize', fontSize)
set(gca, 'FontSize', fontSize, 'FontWeight', 'b')
xlabel('x','FontSize', fontSize, 'FontWeight', 'b')
ylabel('y','FontSize', fontSize, 'FontWeight', 'b')

print('-djpeg', strcat(reconstructionFile,'.jpg'), '-r360');
save(strcat(reconstructionFile,'.mat'), 'x', 'y', 'Reconstruction') ;

end

