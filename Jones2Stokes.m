function [pdI, pdQ, pdU, pdV] = Jones2Stokes(pdH, pdV)
% function [pdI, pdQ, pdU, pdV] = Jones2Stokes(pdH, pdV)
% bhp 20070510

pdI =    pdH.*conj(pdH) + pdV.*conj(pdV);
pdQ =    pdH.*conj(pdH) - pdV.*conj(pdV);
pdU =    pdH.*conj(pdV) + conj(pdH).*pdV;
pdV = i*(pdH.*conj(pdV) - conj(pdH).*pdV);
