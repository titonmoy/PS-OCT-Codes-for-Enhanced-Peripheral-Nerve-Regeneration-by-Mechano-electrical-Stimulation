function [pdAI, pdAQ, pdAU, pdAV] = AverageStokes(pdI, pdQ, pdU, pdV, nPoints, nAlines)
% function [pdAI, pdAQ, pdAU, pdAV] = AverageStokes(pdI, pdQ, pdU, pdV, nPoints, nAlines)
% bhp 20070510

pdFilter = ones([nPoints, nAlines]) / (nPoints*nAlines);
pdAI = imfilter(pdI, pdFilter, 'replicate');
pdAQ = imfilter(pdQ, pdFilter, 'replicate');
pdAU = imfilter(pdU, pdFilter, 'replicate');
pdAV = imfilter(pdV, pdFilter, 'replicate');
pdFactor = pdAI ./ sqrt(pdAQ.^2+pdAU.^2+pdAV.^2);
pdAQ = pdFactor.*pdAQ;
pdAU = pdFactor.*pdAU;
pdAV = pdFactor.*pdAV;

clear pdFilter pdFactor;
