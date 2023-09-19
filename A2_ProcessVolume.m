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

strSaveDir = strcat(strIMAQDir, '\ProcessedWithStokes\');
mkdir(strSaveDir)


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

nPoints = 3;
nAlines = 5;
pdFilter = ones([nPoints, nAlines]) / (nPoints*nAlines);
pdSettings = [nPoints, nAlines];
bAverageSurface = 0;

nFrames = 200;


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


%% load surface

strFileName = strcat('MatFiles\NerveSurfaceTop\', strSample, '_', strScan, '_NerveSurfaceTop');
load(strFileName)


%%


for i = 1 : 1
    %% running reference

    pdReference = zeros(nLineLength, nNumberLines, 2);
    for j = 1 : nFrames
        nIndex = (i-1) * nFrames + j;
        strFilename = listIMAQ(nIndex).name;
        strFilepath = strcat(strIMAQDir, '\', strFilename);

        cellArrays = readHeader(strFilepath);
        [pdIMAQ, ~] = readData(strFilepath, cellArrays);

        pdReference = pdReference + pdIMAQ;
    end
    
    pdReference = pdReference / nFrames;
    
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
    
    for j = 1 : nFrames
        %% read file
        nIndex = (i-1) * nFrames + j;
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

        pdIntensityParallel = (pdIntensityParallelOdd + pdIntensityParallelEven) / 2;
        
        pdIntensityOdd = pdIntensityParallelOdd + pdIntensityPerpendicularOdd;
        pdIntensityEven = pdIntensityParallelEven + pdIntensityPerpendicularEven;
        pdIntensity = (pdIntensityOdd + pdIntensityEven) / 2;


        %% dopu, cumulative phase retardation, and optic axis
        
        pnSurface = pnSurfaceVolume(j,:);
        [pdOpticAxis, pdPhaseRetardation, pdDOPU] = VectorBirefDOPU(pcdParallelOdd, pcdPerpendicularOdd, ...
            pcdParallelEven, pcdPerpendicularEven,  pdSettings, pnSurface, bAverageSurface);
        
        
        %% save
        
        save(strcat(strSaveDir, strFilename(1:end-4)), 'pdIntensityParallel', 'pdIntensity', ...
            'pdDOPU', 'pdPhaseRetardation', 'pdOpticAxis')
        
        disp(nIndex)
    end
end