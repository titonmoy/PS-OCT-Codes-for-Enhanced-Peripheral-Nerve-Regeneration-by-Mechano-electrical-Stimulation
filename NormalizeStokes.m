function [pdNormalized] = NormalizeStokes(pdS)
% function [pdNormalized] = normalizeStokes(pdS)
% bhp 20070510

pdMag = sqrt(dot(pdS, pdS, 3));
pdZeros = find(pdMag == 0.0);
pdMag(pdZeros) = 1.0;
pdQ = pdS(:, :, 1) ./ pdMag;
pdU = pdS(:, :, 2) ./ pdMag;
pdV = pdS(:, :, 3) ./ pdMag;
pdQ(pdZeros) = 0.0;
pdU(pdZeros) = 0.0;
pdV(pdZeros) = 0.0;
pdNormalized(:, :, 3) = pdV;
pdNormalized(:, :, 2) = pdU;
pdNormalized(:, :, 1) = pdQ;
clear pdMag pdZeros pdQ pdU pdV;
