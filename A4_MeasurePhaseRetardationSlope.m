clc
close all
clear


%% data folders

% location of the data

strSample = 'Rat8';


%% image parameters

nHeight = 512;
nWidth = 512;
nFrames = 200;
nSections = 5;

nPoints = 1;
nAlines = 5;
pdFilter = ones([nPoints, nAlines]) / (nPoints*nAlines);

nPoints = 1;
nAlines = 5;

nSlopeMeasuringPoints = 40;

nX = 5;
nZ = 3;
pdFilterEnface = ones([nX, nZ]) / (nX*nZ);

nSmoothingKernel = 10;

dIScale1 = 50;
dIScale2 = 90;

dDOPScale1 = 0.5;
dDOPScale2 = 1;

dPRScale1 = 0.0;
dPRScale2 = 0.01;

nSlopeEdgePolyDegree = 4;
nSlopeEdgeOffset = 5;


%% variables

pdIEnface = NaN(nWidth, nFrames, nSections);
pdDOPUEnface = NaN(nWidth, nFrames, nSections);
pdPRSlopeEnface = NaN(nWidth, nFrames, nSections);


%% manually detected slope measurement edges

pnFrames = [1 10:10:200];
pnY1Sections = zeros(nSections, length(pnFrames));
pnY2Sections = zeros(nSections, length(pnFrames));


pnY1Sections(1, :) = [114, 113, 110, 122, 107, 83, 77, 88, 51, 59, 65, ...
    69, 86, 86, 94, 82, 47, 50, 60, 73, 59];

pnY2Sections(1, :) = [165, 164, 150, 186, 155, 148, 106, 120, 155, 166, 144, ...
    165, 195, 219, 213, 218, 179, 181, 176, 219, 233];


pnY1Sections(2, :) = [51, 51, 33, 34, 65, 33, 56, 22, 53, 29, 58, ...
    11, 42, 34, 36, 22, 12, 10, 14, 14, 14];

pnY2Sections(2, :) = [198, 198, 183, 199, 180, 174, 181, 134, 180, 165, 156, ...
    218, 149, 184, 136, 116, 102, 111, 107, 107, 107];


pnY1Sections(3, :) = [45, 54, 69, 62, 52, 55, 39, 20, 23, 16, 76, ...
    89, 58, 96, NaN, 42, NaN, 67, 52, 61, 108];

pnY2Sections(3, :) = [130, 134, 127, 116, 110, 121, 90, 111, 121, 94, 192, ...
    184, 158, 188, NaN, 145, NaN, 137, 123, 156, 208];


pnY1Sections(4, :) = [51, 65, 59, 89, 83, 112, 101, 103, 79, 65, 62, ...
    59, 41, 46, 41, 43, 41, 44, 45, 49, 40];

pnY2Sections(4, :) = [126, 142, 165, 194, 173, 190, 187, 176, 157, 153, 130, ...
    128, 146, 142,124, 102, 129, 130, 132, 140, 115];


pnY1Sections(5, :) = [61, 76, 47, 44, 66, 71, 68, 96, 77, 103, 64, ...
    67, 78, 95, 43, 66, NaN, 214, 213, 206, 206];

pnY2Sections(5, :) = [132, 147, 112, 126, 125, 124, 126, 180, 173, 165, 167, ...
    149, 188, 156, 140, 155, NaN, 329, 262, 278, 278];



%%

for l = 1 : nSections
    strScan = strcat('Section', num2str(l));

    strDataFolder = strcat(strSample, '\', strScan);

    strMatDir = strcat('..\', strDataFolder, '\ProcessedWithStokesSurfaceHeightAdjusted\');
    listMat = dir(strcat(strMatDir, '\*.mat'));

    strFileName = strcat('MatFiles\NerveSurfaceTop\', strSample, '_', strScan, '_NerveSurfaceTop');
    load(strFileName)

    pnFrameCopy = pnFrames;
    pnY1 = pnY1Sections(l, :);
    pnY2 = pnY2Sections(l, :);

    TF = isnan(pnY1);
    pnFrameCopy(TF) = [];
    pnY1(TF) = [];
    pnY2(TF) = [];

    pnY1Fit = round(interp1(pnFrameCopy, pnY1, 1:nFrames));
    pnY2Fit = round(interp1(pnFrameCopy, pnY2, 1:nFrames));

    for i = 1 : 1
        for j = 1 : nFrames
            nIndex = (i-1) * nFrames + j;
            strFilename = listMat(nIndex).name;
            strFilepath = strcat(strMatDir, '\', strFilename);
            load(strFilepath)

            pdPhaseRetardation = imfilter(pdPhaseRetardation, pdFilter, 'replicate');

            pnSurface = pnSurfaceVolume(j,:);
            pnValidLines = find(pnSurface>1);

            nY1 = pnY1Fit(j) + nSlopeEdgeOffset;
            nY2 = pnY2Fit(j) - nSlopeEdgeOffset;

            for k = pnValidLines
                nLineLength = nHeight - pnSurface(k) + 1;

                pdILine = pdIntensity(:, k);
                pdDOPLine = pdDOPU(:, k);
                pdPRLine = pdPhaseRetardation(:, k);

                pdIEnface(k, j, l) = mean(pdILine(nY1:nY2));
                pdDOPUEnface(k, j, l) = mean(pdDOPLine(nY1:nY2));
                
                p = polyfit(nY1:nY2, pdPRLine(nY1:nY2), 1);
                pdPRSlopeEnface(k, j, l) = p(1);
            end
        end
    end
end


%%

strFileName = strcat('MatFiles\EnFace\', strSample, '_EnFace');
save(strFileName, 'pdIEnface', 'pdDOPUEnface', 'pdPRSlopeEnface')