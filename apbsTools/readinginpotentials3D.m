function [outim3D, thickness, dxnew, dynew, dthick] = readinginpotentials3D(filename)
%readinginpotentials3D: parses the output of APBS and makes a volume from it
% SYNOPSIS:
% [outim3D, thickness, dxnew, dynew, dthick] = readinginpotentials3D(filename)
%
% PARAMETERS:
%   filename: name of APBS file
%
% OUTPUT:
%    outim3D: output volume
%  thickness: thickness of the volume
%      dxnew: voxel size in x direction
%      dynew: voxel size in y direction
%     dthick: voxel size in z direction

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

fid = fopen(filename,'r');% Open text file

%Read Introduction Lines
InputText = textscan(fid,'%s',4,'delimiter','\n'); % Read strings delimited by
%a carriage return
Intro = InputText{1};
disp(Intro);
%Read Each Block
%For each block, we read a header, a table name, column headers for the data, then the data itself.
Block = 1;   % Initialize block index
    sprintf('Block: %s', num2str(Block)); % Display block number
    %InputText=textscan(fid,'%s',1,'delimiter','\n'); % Read header line
    InputText=textscan(fid,'object 1 class gridpositions counts %f %f %f');
    xgrid=InputText{1};
    ygrid=InputText{2};
    zgrid=InputText{3};
Block = 2;
    sprintf('Block: %s', num2str(Block)); % Display block numberot_
    InputText=textscan(fid,'rigin %f %f %f');
    xpos=InputText{1};
    ypos=InputText{2};
    zpos=InputText{3};
    origin=[ zpos ypos xpos];
Block = 3;
    sprintf('Block: %s', num2str(Block)); % Display block number
    InputText=textscan(fid,'delta %f %f %f');
    xdis=InputText{1};
    ydis=InputText{2};
    zdis=InputText{3};
    spacing = [xdis(1) ydis(2)  zdis(3)];
Block = 4;
    sprintf('Block: %s', num2str(Block)); % Display block number
    InputText=textscan(fid,'%s',1,'delimiter','\n'); % Read header line
    HeaderLines{Block,1}=InputText{1};
    disp(HeaderLines{Block});
Block = 5;
    sprintf('Block: %s', num2str(Block)); % Display block number
    InputText=textscan(fid,'%s',1,'delimiter','\n'); % Read header line
    HeaderLines{Block,1}=InputText{1};
    disp(HeaderLines{Block});
% Table of data
Block = 6;
    NumCols=3;
    FormatString=repmat('%f',1,NumCols); % Create format string based on parameter
    InputText=textscan(fid,FormatString,'delimiter',','); % Read data block
    Data{1,1}=cell2mat(InputText); % Convert to numerical array from cell
    [NumRows,NumCols]=size(Data{1});  % Size of table
    disp(cellstr([xlate('Table data size: ') num2str(NumRows) ' x ' num2str(NumCols)]));
    disp(' '); % New line
    eob=textscan(fid,'%s',5,'delimiter','\n');  % Read and discard EOB marker ('EOF' in this case)
%Close the Text File
fclose(fid);

g = Data{1};%matrix with data
nx= xgrid;%=193; % read from .dx file
ny= ygrid;%=289;
nz= zgrid;%=417;
vel=size(g,1)*size(g,2);
kkk=mod(xgrid*ygrid*zgrid,3);
gvec=reshape(g',1,vel);

if kkk==0
 gcut=gvec(1:end);
elseif kkk==1;
 gcut=gvec(1:end-2);
else 
 gcut=gvec(1:end-1);
end

gr = reshape(gcut,zgrid, ygrid, xgrid); %z is the fastest one
grper = gr;
dthick = spacing(1);
dxnew  = spacing(3);
dynew  = spacing(2);
minn   = nx;
thickness = dthick*minn*1e-10;
outim3D = dip_image(grper);



