clc
close all
clear


%% data folders

% location of the data

strSample = 'Rat8';

strFileName = strcat('MatFiles\EnFace\', strSample, '_EnFace');
load(strFileName)


%% parameters

nHr = 39;
nWr = 15;

nHeight = 512;
nWidth = 512;
nFrames = 200;
nSections = 5;

nX = 9;
nZ = 5;
pdFilterEnface = ones([nX, nZ]) / (nX*nZ);

nROIdistanceY = 5;
nROIdistanceX = 2;


%% filtering

pdIEnfaceFiltered = imfilter(pdIEnface, pdFilterEnface);
pdDOPUEnfaceFiltered = imfilter(pdDOPUEnface, pdFilterEnface);

pdPRSlopeEnfaceFiltered = pdPRSlopeEnface;
pdPRSlopeEnfaceFiltered(pdPRSlopeEnface<0) = 0;
pdPRSlopeEnfaceFiltered = imfilter(pdPRSlopeEnfaceFiltered, pdFilterEnface);


%% measurement region

nSectionP = 1; % proxiaml section
nF1P = 1; % proxiaml section starting frame
nF2P = 70; % proxiaml section ending frame

nSectionD = 5; % distal section
nF1D = 165; % distal section staring frame
nF2D = 200; % distal section ending frame



%% variables

pnX1P = []; % proximal tile left edges
pnX2P = []; % proximal tile right edges
pnY1P = []; % proximal tile top edges
pnY2P = []; % proximal tile bottom edges

pnX1D = []; % distal tile left edges
pnX2D = []; % distal tile right edges
pnY1D = []; % distal tile top edges
pnY2D = []; % distal tile bottom edges

pdBirefMP = []; % proximal tile mean of phase retardation slopes
pdBirefSP = []; % proximal tile std of phase retardation slopes

pdBirefMD = []; % distal tile mean of phase retardation slopes
pdBirefSD = []; % distal tile mean of phase retardation slopes


%% proximal

nSection = nSectionP;
nF1 = nF1P;
nF2 = nF2P;

strScan = strcat('Section', num2str(nSection));
strFileName = strcat('MatFiles\NerveSurfaceTop\', strSample, '_', strScan, '_NerveSurfaceTop');
load(strFileName)
pnSurfaceVolume = pnSurfaceVolume';

pdPR = pdPRSlopeEnfaceFiltered(:,:,nSection);

pnMask = zeros(size(pnSurfaceVolume));
pnMask(pnSurfaceVolume>1) = 1;

for j = nF1 : nWr+nROIdistanceX : nF2-nWr
    nX1 = j;
    nX2 = nX1 + nWr - 1;

    pnValidY = find(pnMask(:, j)==1);

    for i = min(pnValidY) : nHr+nROIdistanceY : max(pnValidY)-nHr
        nY1 = i;
        nY2 = nY1 + nHr - 1;

        pnMaskROI = pnMask(nY1:nY2, nX1:nX2);
        if sum(pnMaskROI(:)) < nHr*nWr
            continue
        end

        pnX1P = [pnX1P nX1];
        pnX2P = [pnX2P nX2];

        pnY1P = [pnY1P nY1];
        pnY2P = [pnY2P nY2];

        pdPRroi = pdPR(nY1:nY2, nX1:nX2);
        pdBirefMP = [pdBirefMP mean(pdPRroi(:), 'omitnan')];
        pdBirefSP = [pdBirefSP std(pdPRroi(:), 'omitnan')];

        nCount = nCount+1;
    end
end


%% distal

nSection = nSectionD;
nF1 = nF1D;
nF2 = nF2D;

strScan = strcat('Section', num2str(nSection));
strFileName = strcat('MatFiles\NerveSurfaceTop\', strSample, '_', strScan, '_NerveSurfaceTop');
load(strFileName)
pnSurfaceVolume = pnSurfaceVolume';

pdPR = pdPRSlopeEnfaceFiltered(:,:,nSection);

pnMask = zeros(size(pnSurfaceVolume));
pnMask(pnSurfaceVolume>1) = 1;

for j = nF1 : nWr+nROIdistanceX : nF2-nWr
    nX1 = j;
    nX2 = nX1 + nWr - 1;

    pnValidY = find(pnMask(:, j)==1);

    for i = min(pnValidY) : nHr+nROIdistanceY : max(pnValidY)-nHr
        nY1 = i;
        nY2 = nY1 + nHr - 1;        

        pnMaskROI = pnMask(nY1:nY2, nX1:nX2);
        if sum(pnMaskROI(:)) < nHr*nWr
            continue
        end

        pnX1D = [pnX1D nX1];
        pnX2D = [pnX2D nX2];

        pnY1D = [pnY1D nY1];
        pnY2D = [pnY2D nY2];

        pdPRroi = pdPR(nY1:nY2, nX1:nX2);
        pdBirefMD = [pdBirefMD mean(pdPRroi(:), 'omitnan')];
        pdBirefSD = [pdBirefSD std(pdPRroi(:), 'omitnan')];
    end
end


%% save

strFileName = strcat('MatFiles\EnFace\', strSample, '_BirefProximalDistal');
save(strFileName, 'nSectionP', 'nSectionD', 'nF1P', 'nF2P', 'nF1D', 'nF2D', 'pnX1P', 'pnX2P', 'pnY1P', 'pnY2P', ...
    'pnX1D', 'pnX2D', 'pnY1D', 'pnY2D', 'pdBirefMP', 'pdBirefSP', 'pdBirefMD', 'pdBirefSD')