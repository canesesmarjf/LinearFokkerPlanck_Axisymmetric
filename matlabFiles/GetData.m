% The objective of this script is to extract the binary data from the files
% and create variables for post-processing
clear all
close all

% Define location of data:
% =========================================================================
folderName = '2021_01_27a';
homeAddress = cd;
targetAddress = [homeAddress,'\',folderName];

% Extract simulation conditions:
% =========================================================================
fileName = [targetAddress,'\data.out'];
metadata = GetMetadata(fileName,1);
dt  = metadata.DT;
Te0 = metadata.TE0;
Ti0 = metadata.TI0;
Ne0 = metadata.NE0;
try
    zmin = metadata.ZDUMP;
    zmax = metadata.ZTARGET;
catch
    zmin = metadata.ZMIN;
    zmax = metadata.ZMAX;
end

cd(targetAddress)
% Extract Binary data:
% =========================================================================
% Create file IDs:
% =========================================================================
fid{1}  = fopen('zp.out','r');
fid{2}  = fopen('kep.out','r');
fid{3}  = fopen('xip.out','r');
fid{4}  = fopen('tp.out','r');
fid{5}  = fopen('pcount1.out','r');
fid{6}  = fopen('pcount2.out','r');
fid{7}  = fopen('pcount3.out','r');
fid{8}  = fopen('pcount4.out','r');
fid{9}  = fopen('ecount1.out','r');
fid{10} = fopen('ecount2.out','r');
fid{11} = fopen('ecount3.out','r');
fid{12} = fopen('ecount4.out','r');
cd(homeAddress)

% Read binary data:
% =========================================================================
f1 = ExtractBinaryData(fid(1:4),8);
f2 = ExtractBinaryData(fid(5:12),8);

% Assign data to variables:
% =========================================================================
zp        = f1{1};
kep       = f1{2};
xip       = f1{3};
tp        = f1{4};
pcount1_0 = f2{1};
pcount2_0 = f2{2};
pcount3_0 = f2{3};
pcount4_0 = f2{4};
ecount1_0 = f2{5};
ecount2_0 = f2{6};
ecount3_0 = f2{7};
ecount4_0 = f2{8};
tc        = 0:dt:((numel(pcount1_0)-1)*dt);
clearvars f1 f2 fid

% Magnetic field data
% =========================================================================
BfieldAddress = dir([targetAddress,'\*.txt']);
f = load([BfieldAddress.folder,'\',BfieldAddress.name]);
zb = f(:,1);
b  = f(:,2);
bmax = max(b);

% Plot magnetic field data:
figure('color','w')
plot(zb,b)
grid on
box on
xlabel('z [m]')

%% Define functions:
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