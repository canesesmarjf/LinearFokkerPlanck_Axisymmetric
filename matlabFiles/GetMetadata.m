function [metadata] = GetMetadata(fileName,plotType)

% Assign file identifier:
% =========================================================================
fileID = fopen(fileName,'r');

% Scan text:
% =========================================================================
delimiter = {'='};
formatSpec = '%q%[^\n\r]';
dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter);

% Create "metadata" structure:
% =========================================================================
for ii = 2:(numel(dataArray{1})-1)
    dum1 = str2double(dataArray{2}{ii}(1:end-1));
    if isnan(dum1)
       metadata.(dataArray{1}{ii}) = dataArray{2}{ii}(1:end-1);
    else
       metadata.(dataArray{1}{ii}) = dum1;
    end
end

fclose(fileID);

switch plotType
    case 1
        clc
        metadata
    case 2
        clc
        T = struct2table(metadata,'AsArray',1)
end

end

