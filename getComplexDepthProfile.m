function pcdDepthProfile = getComplexDepthProfile(pcdIMAQ, pdMask)

[~, nNumberLines] = size(pcdIMAQ);

pcdIMAQ = pcdIMAQ .* repmat(pdMask, [1 nNumberLines]);
pcdDepthProfile = fft(pcdIMAQ);

end