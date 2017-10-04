%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTERED BACKPROJECTION ALGORITHM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FBP(dataFileName, reconstructionFile, x, y)

fontSize = 18;

load(dataFileName);

% --- Smallest power of 2 greater than numDetectors
n1 = ceil(log(numDetectors) / log(2)); 
N = 2^n1;                  
% N = 2^(n1 + 1);

% --- Inter-detector distance
dt = (max(detectorPositions) - min(detectorPositions)) / (numDetectors - 1);  

% --- Maximum range value after zero padding
tmax = dt * N;             
% --- Sampling step in the frequency domain
dW = 1 / tmax;         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZERO-PADDED FFTs OF THE PROJECTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filteredProjectionSpectra = fft(sinogram.', N);        % --- Padded FFT is performed for every column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTERING IN THE SPECTRAL DOMAIN %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = [0 : (N / 2) (N / 2 - 1) : -1 : 1].' * dW;
filteredProjectionSpectra = w .* filteredProjectionSpectra;

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
title('Filtered BackProjection', 'FontSize', fontSize)
set(gca, 'FontSize', fontSize, 'FontWeight', 'b')
xlabel('x','FontSize', fontSize, 'FontWeight', 'b')
ylabel('y','FontSize', fontSize, 'FontWeight', 'b')

print('-djpeg', strcat(reconstructionFile,'.jpg'), '-r360');
save(strcat(reconstructionFile,'.mat'), 'x', 'y', 'Reconstruction') ;

end

