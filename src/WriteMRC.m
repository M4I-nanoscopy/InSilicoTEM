function WriteMRC(map,pixA,filename,mode,mz)
% function WriteMRC(map,pixA,filename,mode,mz)
%  (To get a file selector, give no filename or give [] as the filename.)
% Write out a 2D image or a 3D volume as an MRC map file, for example for viewing in
% Chimera.  'map' is the 3D array, pixA is the voxel size in angstroms.
% mz is the number of slices, e.g. when writing a stack of 2D images.
% Mode sets the data type:
% 0: uint8
% 1: int16
% 2: single (default)
% 6: uint16
% 32: packed bytes (uint4)  **New**
% fixed for 2D output as well.  fs 28 Aug 07
% % Test program: create an ellipsoid and store it as a map file.
%   [x y z]=ndgrid(-32:31);
%   m=(x.^2+(y*1.5).^2+z.^2)>20^2;
%   WriteMRC(m,1,'maptest.mrc',2,64); % mz is the unit cell size in z.

if nargin<3 || numel(filename)<1
    fprintf('Getting the MRC filename to write: ');
    [name,path]=uiputfile('*.mrc;*.mrcs','Writing an MRC file');
    if isnumeric(name)
        warning('File not written.');
        return
    else
        disp([path name]);
        filename=[path name];
    end;
end;
if nargin<4
    mode=2;
end;
if nargin<5
    mz=1;  % images
end;

sizes=[size(map,1) size(map,2) size(map,3)];
org=-floor(sizes/2);  % Default origin
org=[0 0 0];
switch mode
    case 0
        string='int8';
    case 1
        string='int16';
    case 2
        string='float32';
    case 6
        string='uint16';
        
    case 32  % Packed uint4 values
        string='uint8';
        map=uint8(map);
        tempMap=min(map,15);
        %         The following error checking takes about 30% of total execution time.
        frac=sum(map(:)>tempMap(:))/numel(map);
        if frac>1e-3
            error(['4-bit packed bytes: fraction out of range is ' num2str(frac)]);
        elseif frac>1e-5
            warning(['4-bit packed bytes: fraction out of range is ' num2str(frac)]);
        end;
        %           Make sure we have an even number in each row.
        if rem(sizes(1),2)==1
            tempMap(end+1,:,:)=0; % Force an even number.
        end;
        tempMap=reshape(tempMap,2,numel(tempMap)/2);
        map=uint8(tempMap(1,:)+16*tempMap(2,:));
    otherwise
        error(['Invalid data mode: ' num2str(mode)]);
end;
f=WriteMRCHeader(map,pixA,filename,sizes,org,mode,mz);
count2=fwrite(f,map,string);
fclose(f);

return;
