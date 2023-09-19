clc
close all
clear


%%  data folders

% location of the data

strSample = 'Rat8';
strScan = 'Section1';

strDataFolder = strcat(strSample, '\', strScan);

strMatDir = strcat('..\', strDataFolder, '\ProcessedWithStokes\');
listMat = dir(strcat(strMatDir, '\*.mat'));

strSaveDir = strcat('..\', strDataFolder, '\ProcessedWithStokesSurfaceHeightAdjusted\');
mkdir(strSaveDir)

strFileName = strcat('MatFiles\NerveSurfaceTop\', strSample, '_', strScan, '_NerveSurfaceTop');
load(strFileName)


%% image parameters

nHeight = 512;
nWidth = 512;
nFrames = 200;


%%

pdIntensityParallelNaN = NaN(nHeight, nWidth);
pdIntensityNaN = NaN(nHeight, nWidth);
pdDOPUNaN = NaN(nHeight, nWidth);
pdPhaseRetardationNaN = NaN(nHeight, nWidth);
pdOpticAxisNaN = NaN(nHeight, nWidth, 3);

for i = 1 : 1    
    for j = 1 : nFrames
        nIndex = (i-1) * nFrames + j;
        strFilename = listMat(nIndex).name;
        strFilepath = strcat(strMatDir, '\', strFilename);
        load(strFilepath)

        pnSurface = pnSurfaceVolume(j,:);

        pdIntensityParallelNaN = pdIntensityParallelNaN * NaN;
        pdIntensityNaN = pdIntensityNaN * NaN;
        pdDOPUNaN = pdDOPUNaN * NaN;
        pdPhaseRetardationNaN = pdPhaseRetardationNaN * NaN;
        pdOpticAxisNaN = pdOpticAxisNaN * NaN;

        for k = 1 : nWidth
            nHeightAdjusted = nHeight - pnSurface(k) + 1;
            pdIntensityParallelNaN(1:nHeightAdjusted, k) = pdIntensityParallel(pnSurface(k) : nHeight, k);
            pdIntensityNaN(1:nHeightAdjusted, k) = pdIntensity(pnSurface(k) : nHeight, k);
            pdDOPUNaN(1:nHeightAdjusted, k) = pdDOPU(pnSurface(k) : nHeight, k);
            pdPhaseRetardationNaN(1:nHeightAdjusted, k) = pdPhaseRetardation(pnSurface(k) : nHeight, k);
            pdOpticAxisNaN(1:nHeightAdjusted, k, :) = pdOpticAxis(pnSurface(k) : nHeight, k, :);
        end

        pdIntensityParallel = pdIntensityParallelNaN;
        pdIntensity = pdIntensityNaN;
        pdDOPU = pdDOPUNaN;
        pdPhaseRetardation = pdPhaseRetardationNaN;
        pdOpticAxis = pdOpticAxisNaN;

        save(strcat(strSaveDir, strFilename(1:end-4)), 'pdIntensityParallel', 'pdIntensity', ...
            'pdDOPU', 'pdPhaseRetardation', 'pdOpticAxis')
        
        disp(nIndex)
    end
end