function [f] = ExtractBinaryData(fid,floatSize)
% EXTRACTBINARYDATA:
% =========================================================================
% This function is made specifically for the MPEX_LinearFokkerPlanck
% Fortran code developed in 2019
% Extracts the binary data contained in fID given the size of the floating
% point number it represents in Bytes, where fid is the file identifier
% produced using fID = fopen(fileName,[option1],[option2],...)
% INPUTS:
% =========================================================================
% fid: Cell of N file identifier
% floatSize: (Scalar) Number of Bytes of the floating point numbers
% represented by in the binary file associated with the file identifiers
% OUTPUTS:
% =========================================================================
% f: Cell with N elements. Each element is an array with the dimensions of
% the data it represents
% =========================================================================

% Determine size of arrays:
% -------------------------------------------------------------------------
for ii = 1:numel(fid)
    nSize(ii) = (fread(fid{ii},1,'uint32'))/floatSize;
    frewind(fid{ii})
end
nSizeMin = min(nSize);
nSizeMax = max(nSize);
for ii = 1:numel(fid)
    arraySize{ii} = [nSize(ii)/nSizeMin,nSizeMin];
end

% Read the binary files:
% -------------------------------------------------------------------------
for ii = 1:numel(fid)
    [~] = fread(fid{ii},1,'uint32');  
    f{ii} = fread(fid{ii},arraySize{ii},'real*8');
    fclose(fid{ii});
end
end

