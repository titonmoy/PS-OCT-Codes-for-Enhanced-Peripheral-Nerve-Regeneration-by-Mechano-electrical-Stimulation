clc
close all
clear


%% data folders

% location of the data

strSample = 'Rat8';
strScan = 'Section1';

strDataFolder = strcat(strSample, '\', strScan);

strIMAQDir = strcat('..\', strDataFolder);
listIMAQ = dir(strcat(strIMAQDir, '\*.dat'));

strRefDir = strcat(strIMAQDir, '\Ref');
listRef = dir(strcat(strRefDir, '\*.dat'));


%% parameters

% system and scan specific parameters

stParameter = readtable('..\parameter.txt');
nParameterIndex = find(strcmp(strDataFolder, stParameter.Folder));

nSample = stParameter.Line(nParameterIndex);
nShift = stParameter.LineShift(nParameterIndex);
nSpectrumShift = stParameter.SpectrumShift(nParameterIndex);

nNumberCameras = 2;
nNumberPolarizationStates = 2;
nLineSamplesPerState = 2;
nNumberCalibrationLines = nNumberCameras * nNumberPolarizationStates;

nNumberOCTLinesPerPolarizationState = 512;
nLineLength = 1024;
nNumberLines = nLineSamplesPerState * nNumberPolarizationStates * nNumberOCTLinesPerPolarizationState;

nInterval = nLineSamplesPerState * nNumberPolarizationStates;

nFrames = 200;

nHeight = 512;
nWidth = 512;

nPoints = 3;
nAlines = 5;
pdFilter = ones([nPoints, nAlines]) / (nPoints*nAlines);

dIScale1 = 50;
dIScale2 = 90;

dDOPScale1 = 0;
dDOPScale2 = 1;


%% load calibration

% system specific calibration file

load('calibration.mat')
pdKParallelOdd = pdK(:, 1);
pdKParallelEven = pdK(:, 2);
pdKPerpendicularOdd = pdK(:, 3);
pdKPerpendicularEven = pdK(:, 4);
clear pdK

pnIndex = zeros(nLineLength, nNumberCalibrationLines);
pnIndexParallelOdd = pnIndex(1, :)';
pnIndexParallelEven = pnIndex(2, :)';
pnIndexPerpendicularOdd = pnIndex(3, :)';
pnIndexPerpendicularEven = pnIndex(4, :)';
clear pnIndex

bMATLABInterp = true;


%% read dispersion

% system specific dispersion file

load('dispersion_window.mat')


%% calculate mask

% window generation

nLeft = 100;
nRight = 924;
nRound = 16;

pdMask = calculateMask(nLineLength, nLeft, nRight, nRound);


%% pixel correction parameters

% system specific parameter

nPixelLeft = 414;
nPixelRight = 419;


%% read reference

filename = listRef(1).name;
strRef = strcat(strRefDir, '\', filename);

cellArrays = readHeader(strRef);
[pdReference, ~] = readData(strRef, cellArrays);

pdReferenceParallel = pdReference(:, :, 1);
pdReferencePerpendicular = flipud(pdReference(:, :, 2));

pdReferenceParallelOdd = pdReferenceParallel(:, nSample : nInterval : end);
pdReferenceParallelEven = pdReferenceParallel(:, nSample+nLineSamplesPerState : nInterval : end);

pdReferencePerpendicularOdd = pdReferencePerpendicular(:, nSample : nInterval : end);
pdReferencePerpendicularEven = pdReferencePerpendicular(:, nSample+nLineSamplesPerState : nInterval : end);

pdReferenceParallelOddMean = mean(pdReferenceParallelOdd, 2);
pdReferenceParallelEvenMean = mean(pdReferenceParallelEven, 2);

pdReferencePerpendicularOddMean = mean(pdReferencePerpendicularOdd, 2);
pdReferencePerpendicularEvenMean = mean(pdReferencePerpendicularEven, 2);

% system specific correction

dSlope = (pdReferencePerpendicularOddMean(nPixelRight) - pdReferencePerpendicularOddMean(nPixelLeft)) / (nPixelRight - nPixelLeft);
dOffset = pdReferencePerpendicularOddMean(nPixelLeft) - dSlope * nPixelLeft;
pdReferencePerpendicularOddMean(nPixelLeft+1:nPixelRight-1) = dSlope * (nPixelLeft+1:nPixelRight-1) + dOffset;

dSlope = (pdReferencePerpendicularEvenMean(nPixelRight) - pdReferencePerpendicularEvenMean(nPixelLeft)) / (nPixelRight - nPixelLeft);
dOffset = pdReferencePerpendicularEvenMean(nPixelLeft) - dSlope * nPixelLeft;
pdReferencePerpendicularEvenMean(nPixelLeft+1:nPixelRight-1) = dSlope * (nPixelLeft+1:nPixelRight-1) + dOffset;


%% detect surface

pdIntensitydBVolume = zeros(nHeight, nWidth, nFrames);
pdDOPUVolume = zeros(nHeight, nWidth, nFrames);

for i = 1 : nFrames
    %% read file
    nIndex = i;
    strFilename = listIMAQ(nIndex).name;
    strFilepath = strcat(strIMAQDir, '\', strFilename);

    cellArrays = readHeader(strFilepath);
    [pdIMAQ, ~] = readData(strFilepath, cellArrays);


    %% IMAQ

    pdIMAQParallel = pdIMAQ(:, :, 1);
    pdIMAQPerpendicular = flipud(pdIMAQ(:, :, 2));

    pdIMAQParallelOdd = pdIMAQParallel(:, nSample : nInterval : end);
    pdIMAQParallelEven = pdIMAQParallel(:, nSample+nLineSamplesPerState : nInterval : end);

    pdIMAQPerpendicularOdd = pdIMAQPerpendicular(:, nSample : nInterval : end);
    pdIMAQPerpendicularEven = pdIMAQPerpendicular(:, nSample+nLineSamplesPerState : nInterval : end);


    % system specific correction

    dSlope = (pdIMAQPerpendicularOdd(nPixelRight, :) - pdIMAQPerpendicularOdd(nPixelLeft, :)) / (nPixelRight - nPixelLeft);
    dOffset = pdIMAQPerpendicularOdd(nPixelLeft, :) - dSlope * nPixelLeft;
    pdIMAQPerpendicularOdd(nPixelLeft+1:nPixelRight-1, :) = (nPixelLeft+1:nPixelRight-1)' * dSlope + dOffset;

    dSlope = (pdIMAQPerpendicularEven(nPixelRight, :) - pdIMAQPerpendicularEven(nPixelLeft, :)) / (nPixelRight - nPixelLeft);
    dOffset = pdIMAQPerpendicularEven(nPixelLeft, :) - dSlope * nPixelLeft;
    pdIMAQPerpendicularEven(nPixelLeft+1:nPixelRight-1, :) = (nPixelLeft+1:nPixelRight-1)' * dSlope + dOffset;


    %% subtract reference

    pdIMAQParallelOdd = pdIMAQParallelOdd - repmat(pdReferenceParallelOddMean, [1 nNumberOCTLinesPerPolarizationState]);
    pdIMAQParallelEven = pdIMAQParallelEven - repmat(pdReferenceParallelEvenMean, [1 nNumberOCTLinesPerPolarizationState]);

    pdIMAQPerpendicularOdd = pdIMAQPerpendicularOdd - repmat(pdReferencePerpendicularOddMean, [1 nNumberOCTLinesPerPolarizationState]);
    pdIMAQPerpendicularEven = pdIMAQPerpendicularEven - repmat(pdReferencePerpendicularEvenMean, [1 nNumberOCTLinesPerPolarizationState]);


    %% apply calibration

    pdIMAQParallelOdd = applyCalibration(pdIMAQParallelOdd, pdKParallelOdd, pnIndexParallelOdd, bMATLABInterp);

    pdIMAQParallelEven = applyCalibration(pdIMAQParallelEven, pdKParallelEven, pnIndexParallelEven, bMATLABInterp);

    pdIMAQPerpendicularOdd = applyCalibration(pdIMAQPerpendicularOdd, pdKPerpendicularOdd, pnIndexPerpendicularOdd, bMATLABInterp);

    pdIMAQPerpendicularEven = applyCalibration(pdIMAQPerpendicularEven, pdKPerpendicularEven, pnIndexPerpendicularEven, bMATLABInterp);


    %% apply dispersion

    pcdIMAQParallelOdd = applyDispersion(pdIMAQParallelOdd, pdDispersionReal, pdDispersionImag);
    pcdIMAQParallelEven = applyDispersion(pdIMAQParallelEven, pdDispersionReal, pdDispersionImag);
    pcdIMAQPerpendicularOdd = applyDispersion(pdIMAQPerpendicularOdd, pdDispersionReal, pdDispersionImag);
    pcdIMAQPerpendicularEven = applyDispersion(pdIMAQPerpendicularEven, pdDispersionReal, pdDispersionImag);


    %% spectrum shift; system specific correction

    pcdIMAQParallelOdd = circshift(pcdIMAQParallelOdd, -ceil(nSpectrumShift/2), 1);
    pcdIMAQParallelEven = circshift(pcdIMAQParallelEven, -ceil(nSpectrumShift/2), 1);
    pcdIMAQPerpendicularOdd = circshift(pcdIMAQPerpendicularOdd, floor(nSpectrumShift/2), 1);
    pcdIMAQPerpendicularEven = circshift(pcdIMAQPerpendicularEven, floor(nSpectrumShift/2), 1);


    %% shifts for camera sync; system specific correction

    pcdIMAQPerpendicularOdd = circshift(pcdIMAQPerpendicularOdd, nShift, 2);
    pcdIMAQPerpendicularEven = circshift(pcdIMAQPerpendicularEven, nShift, 2);


    %% get complex depth profile

    pcdParallelOdd = getComplexDepthProfile(pcdIMAQParallelOdd, pdMask);
    pcdParallelEven = getComplexDepthProfile(pcdIMAQParallelEven, pdMask);
    pcdPerpendicularOdd = getComplexDepthProfile(pcdIMAQPerpendicularOdd, pdMask);
    pcdPerpendicularEven = getComplexDepthProfile(pcdIMAQPerpendicularEven, pdMask);

    pcdParallelOdd(nLineLength/2+1 : end, :) = [];
    pcdParallelEven(nLineLength/2+1 : end, :) = [];
    pcdPerpendicularOdd(nLineLength/2+1 : end, :) = [];
    pcdPerpendicularEven(nLineLength/2+1 : end, :) = [];


    %% intensity image

    pdIntensityParallelOdd = abs(pcdParallelOdd).^2;
    pdIntensityParallelEven = abs(pcdParallelEven).^2;
    pdIntensityPerpendicularOdd = abs(pcdPerpendicularOdd).^2;
    pdIntensityPerpendicularEven = abs(pcdPerpendicularEven).^2;

    pdIntensityOdd = pdIntensityParallelOdd + pdIntensityPerpendicularOdd;
    pdIntensityEven = pdIntensityParallelEven + pdIntensityPerpendicularEven;
    pdIntensity = (pdIntensityOdd + pdIntensityEven) / 2;
    pdIntensitydB = 10 * log10(pdIntensity);


    %% DOPU

    [pdIE, pdQE, pdUE, pdVE] = Jones2Stokes(pcdParallelEven, pcdPerpendicularEven);
    [pdIO, pdQO, pdUO, pdVO] = Jones2Stokes(pcdParallelOdd, pcdPerpendicularOdd);

    pdIE = imfilter(pdIE, pdFilter, 'replicate');
    pdQE = imfilter(pdQE, pdFilter, 'replicate');
    pdUE = imfilter(pdUE, pdFilter, 'replicate');
    pdVE = imfilter(pdVE, pdFilter, 'replicate');

    pdIO = imfilter(pdIO, pdFilter, 'replicate');
    pdQO = imfilter(pdQO, pdFilter, 'replicate');
    pdUO = imfilter(pdUO, pdFilter, 'replicate');
    pdVO = imfilter(pdVO, pdFilter, 'replicate');

    pdIENorm = pdIE ./ sqrt(pdQE.^2 + pdUE.^2 + pdVE.^2);
    pdIONorm = pdIO ./ sqrt(pdQO.^2 + pdUO.^2 + pdVO.^2);

    pdDOPU = 0.5 * (1./pdIENorm + 1./pdIONorm);


    %% 

    pdIntensitydBVolume(:,:,i) = pdIntensitydB;
    pdDOPUVolume(:,:,i) = pdDOPU;

    disp(i)
end


%% 

pnFrames =  [1 5:5:200];
nNumberOfFrames = length(pnFrames);
pnFrameSurface = ones(nNumberOfFrames, nWidth);


%% manually draw surface

figure
set(gcf, 'Position', get(0, 'Screensize'))

for i = 1 : nNumberOfFrames
    nIndex = pnFrames(i);
    pdI = pdIntensitydBVolume(:,:,nIndex);
    pdDOP = pdDOPUVolume(:,:,nIndex);

    pdI(pdI<dIScale1) = dIScale1;
    pdI(pdI>dIScale2) = dIScale2;
    pdI = (pdI - dIScale1) / (dIScale2 - dIScale1);

    pdDOP(pdDOP<dDOPScale1) = dDOPScale1;
    pdDOP(pdDOP>dDOPScale2) = dDOPScale2;
    pdDOP = (pdDOP - dDOPScale1) / (dDOPScale2 - dDOPScale1);

    pdCompositeHSV = zeros(nHeight, nWidth, 3);
    pdCompositeHSV(:,:,1) = pdDOP;
    pdCompositeHSV(:,:,2) = 1;
    pdCompositeHSV(:,:,3) = pdI;
    pnCompositeRGB = hsv2rgb(pdCompositeHSV);

    clf
    imagesc(pnCompositeRGB)
    ax = gca;
    title({[strSample, ' ', strScan]; ['Frame: ', num2str(nIndex)]})
    hold on
    h = images.roi.Freehand(ax, 'Color', 'r');
    h.draw
    pdEdgePosition = h.Position;
    plot(ax, pdEdgePosition(:,1), pdEdgePosition(:,2), '-m', 'LineWidth', 3)

    pdEdgePositionRound = round(pdEdgePosition(:,1));
    [~, idx, ~] = unique(pdEdgePositionRound(:,1));
    pdEdgePosition = pdEdgePosition(idx, :);
    pnLines = ceil(pdEdgePosition(1,1)) : floor(pdEdgePosition(end,1));
    pnTop = interp1(pdEdgePosition(:,1), pdEdgePosition(:,2), pnLines);
    pnTop = round(pnTop);
    plot(ax, pnLines, pnTop, '-k', 'LineWidth', 3)
    drawnow

    pnFrameSurface(i, pnLines) = pnTop;
end


%% fit surface

pnFrameSurfaceSmooth = zeros(nNumberOfFrames, nWidth);
nSmoothingKernel = 5;
for i = 1 : nNumberOfFrames
    pnFrameSurfaceSmooth(i, :) = smooth(pnFrameSurface(i, :), nSmoothingKernel);
end

pnSurfaceVolume = zeros(nFrames, nWidth);
for j = 1 : nWidth
    pnSurfaceVolume(:, j) = interp1(pnFrames, pnFrameSurfaceSmooth(:, j), 1:nFrames);
end

nSmoothingKernel = 5;
for i = 1 : nFrames
    pnSurfaceVolume(i, :) = smooth(pnSurfaceVolume(i, :), nSmoothingKernel);
end

pnSurfaceVolume = round(pnSurfaceVolume);
pnSurfaceVolume(pnSurfaceVolume<1) = 1;

for i = 1 : nFrames
    pnSurface = pnSurfaceVolume(i,:);
    pnMask = zeros(size(pnSurface));
    pnMask(pnSurface>1) = 1;
    bw = bwareafilt(logical(pnMask), 1);
    pnSurface(~bw) = 1;
    pnSurfaceVolume(i,:) = pnSurface;
end


%% save file

strFileName = strcat('MatFiles\NerveSurfaceTop\', strSample, '_', strScan, '_NerveSurfaceTop');
save(strFileName, 'pnSurfaceVolume')