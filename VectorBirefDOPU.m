function [pdA, pdBiref, pdDOPU] = VectorBirefDOPU(pcdParallelOdd, pcdPerpendicularOdd, pcdParallelEven, pcdPerpendicularEven, pdSettings, pnSurface, bAverageSurface)

% determine Stokes parameters
[pdIE, pdQE, pdUE, pdVE] = Jones2Stokes(pcdParallelEven, pcdPerpendicularEven);
[pdIO, pdQO, pdUO, pdVO] = Jones2Stokes(pcdParallelOdd, pcdPerpendicularOdd);

% average polarization states
nPoints = pdSettings(1);
nAlines = pdSettings(2);
pdFilter = ones([nPoints, nAlines]) / (nPoints*nAlines);

pdIE = imfilter(pdIE, pdFilter, 'replicate');
pdQE = imfilter(pdQE, pdFilter, 'replicate');
pdUE = imfilter(pdUE, pdFilter, 'replicate');
pdVE = imfilter(pdVE, pdFilter, 'replicate');

pdIO = imfilter(pdIO, pdFilter, 'replicate');
pdQO = imfilter(pdQO, pdFilter, 'replicate');
pdUO = imfilter(pdUO, pdFilter, 'replicate');
pdVO = imfilter(pdVO, pdFilter, 'replicate');

% DOPU
pdIENorm = pdIE ./ sqrt(pdQE.^2 + pdUE.^2 + pdVE.^2);
pdIONorm = pdIO ./ sqrt(pdQO.^2 + pdUO.^2 + pdVO.^2);

pdDOPU = 0.5 * (1./pdIENorm + 1./pdIONorm);

pdQE = pdIENorm .* pdQE;
pdUE = pdIENorm .* pdUE;
pdVE = pdIENorm .* pdVE;

pdQO = pdIONorm .* pdQO;
pdUO = pdIONorm .* pdUO;
pdVO = pdIONorm .* pdVO;

% surface states
pdSIE = zeros(1, length(pnSurface));
pdSQE = zeros(1, length(pnSurface));
pdSUE = zeros(1, length(pnSurface));
pdSVE = zeros(1, length(pnSurface));
for nAline = 1 : length(pnSurface);
    pdSIE(1, nAline) = pdIE(pnSurface(nAline), nAline);
    pdSQE(1, nAline) = pdQE(pnSurface(nAline), nAline);
    pdSUE(1, nAline) = pdUE(pnSurface(nAline), nAline);
    pdSVE(1, nAline) = pdVE(pnSurface(nAline), nAline);
end
pdSIO = zeros(1, length(pnSurface));
pdSQO = zeros(1, length(pnSurface));
pdSUO = zeros(1, length(pnSurface));
pdSVO = zeros(1, length(pnSurface));
for nAline = 1 : length(pnSurface);
    pdSIO(1, nAline) = pdIO(pnSurface(nAline), nAline);
    pdSQO(1, nAline) = pdQO(pnSurface(nAline), nAline);
    pdSUO(1, nAline) = pdUO(pnSurface(nAline), nAline);
    pdSVO(1, nAline) = pdVO(pnSurface(nAline), nAline);
end

% average states
if (bAverageSurface);
    pdMQ = mean(pdSQE);  pdMU = mean(pdSUE);  pdMV = mean(pdSVE);
    pdMI = sqrt(pdMQ.*2 + pdMU.^2 + pdMV.^2);
    pdSQE = pdMQ * ones([1, length(pdSIE)]);
    pdSUE = pdMU * ones([1, length(pdSIE)]);
    pdSVE = pdMV * ones([1, length(pdSIE)]);
    clear pdMQ pdMU pdMV pdMI;

    pdMQ = mean(pdSQO);  pdMU = mean(pdSUO);  pdMV = mean(pdSVO);
    pdMI = sqrt(pdMQ.*2 + pdMU.^2 + pdMV.^2);
    pdSQO = pdMQ * ones([1, length(pdSIO)]);
    pdSUO = pdMU * ones([1, length(pdSIO)]);
    pdSVO = pdMV * ones([1, length(pdSIO)]);
    clear pdMQ pdMU pdMV pdMI;
end

% assemble surface Stokes 3-vectors
pdSE(:,:,3) = pdSVE;
pdSE(:,:,2) = pdSUE;
pdSE(:,:,1) = pdSQE;
pdSE = NormalizeStokes(pdSE);
pdSE = repmat(pdSE, [size(pdIE, 1), 1, 1]);
pdSIE = repmat(pdSIE, [size(pdIE, 1), 1]);
pdSO(:,:,3) = pdSVO;
pdSO(:,:,2) = pdSUO;
pdSO(:,:,1) = pdSQO;
pdSO = NormalizeStokes(pdSO);
pdSO = repmat(pdSO, [size(pdIO, 1), 1, 1]);
pdSIO = repmat(pdSIO, [size(pdIO, 1), 1]);
clear pdSQE pdSUE pdSVE pdSQO pdSUO pdSVO;

% assemble Stokes 3-vectors
pdE(:,:,3) = pdVE;
pdE(:,:,2) = pdUE;
pdE(:,:,1) = pdQE;
pdE = NormalizeStokes(pdE);
pdO(:,:,3) = pdVO;
pdO(:,:,2) = pdUO;
pdO(:,:,1) = pdQO;
pdO = NormalizeStokes(pdO);
clear pdQE pdUE pdVE pdQO pdUO pdVO;

% calculate optic axis
pdPE = pdE - pdSE;
pdPO = pdO - pdSO;
pdA = cross(pdPE, pdPO, 3);
pdA = NormalizeStokes(pdA);
clear pdPE pdPO;

% calculate even angles
pdP       = cross(pdA, pdE,  3);
pdSP      = cross(pdA, pdSE, 3);
pdWeightE = pdIE.*dot(pdP, pdP, 3) .* pdSIE.*dot(pdSP, pdSP, 3);
pdP       = NormalizeStokes(pdP);
pdSP      = NormalizeStokes(pdSP);
pdAngleE  = acos(dot(pdP, pdSP, 3));
clear pdP pdSP;

% calculate odd angles
pdP       = cross(pdA, pdO,  3);
pdSP      = cross(pdA, pdSO, 3);
pdWeightO = pdIO.*dot(pdP, pdP, 3) .* pdSIO.*dot(pdSP, pdSP, 3);
pdP       = NormalizeStokes(pdP);
pdSP      = NormalizeStokes(pdSP);
pdAngleO  = acos(dot(pdP, pdSP, 3));
clear pdP pdSP;

% calculate overall angle
clear pdIE pdE pdIO pdO pdSIE pdSE pdSIO pdSO;
pdWeightedAngles = pdWeightE.*pdAngleE + pdWeightO.*pdAngleO;
pdCombinedWeight = pdWeightE           + pdWeightO;
pnZeros          = find(pdCombinedWeight == 0);
pdWeightedAngles(pnZeros) = 0.0;
pdCombinedWeight(pnZeros) = 1.0;
pdBiref          = pdWeightedAngles ./ pdCombinedWeight;

clear pdWeightedAngles pdCombinedWeight;

