function [pdIMAQ, pdDAQ] = readData(strFile, cellArrays)

fp = fopen(strFile, 'r', 'l');
strSysType = cellArrays{1, 3}; 

for nArrayNumber = 2 : size(cellArrays, 1)
    strVar        = cellArrays{nArrayNumber, 1};
    nOffset       = cellArrays{nArrayNumber, 2};
    nNumberLines  = cellArrays{nArrayNumber, 3};
    nNumberPoints = cellArrays{nArrayNumber, 4};
    strDataType   = cellArrays{nArrayNumber, 5};
    fseek(fp, nOffset, 'bof');    
    
    if strcmp(strSysType, 'SD-OCT') % one cameras
        strTest = sprintf('%s = fread(fp, [%d, %d], ''%s'');', strVar, nNumberPoints, nNumberLines, strDataType);
        eval(strTest);
    end 
    
    if strcmp(strSysType, 'PS-SD-OCT') % two cameras
        strTest = sprintf('%s = fread(fp, ''%s'');', strVar, strDataType);
        eval(strTest);
        if strcmp(strVar, 'pdIMAQ')
            pdIMAQ1 = reshape(pdIMAQ(1:2:end), [nNumberPoints, nNumberLines]); % pdIMAQParallel (H)
            pdIMAQ2 = reshape(pdIMAQ(2:2:end), [nNumberPoints, nNumberLines]); % pdIMAQPerpendicular (V)
            pdIMAQ = []; 
            pdIMAQ(:,:,1) = pdIMAQ1; 
            pdIMAQ(:,:,2) = pdIMAQ2;
            clear pdIMAQ1 pdIMAQ2 nIMAQLength; 
        end
    end

    clear strVar nOffset nNumberLines nNumberPoints strDataType;
    clear strTest;
end
fclose(fp);
clear ans fp nArrayNumber;

end