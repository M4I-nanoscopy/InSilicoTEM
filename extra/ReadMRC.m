function [map,s,hdr,extraHeader]=ReadMRC(filename,startSlice, numSlices,test)
% function map=ReadMRC;  % put up a file selector and get a file.  Or,
% function map=ReadMRC(filename);  --all defaults.   Or,
% function [map,s]=ReadMRC(filename,startSlice, numSlices,test)  or,
% function [fh,s]=ReadMRC(filename,startSlice,-1,test)
% Read a 3D map from an .mrc file and return a structure containing various parameters read from the file.
% This function reads 2d and 3d real maps in byte, int16 and float32 modes.
% s.nz gives the total number of slices, and calling [m s]=ReadMRC(name,1,0)
% will give you that value.  startSlice starts at 1.
% The test argument is a boolean.  If true, diagnostic information is
% listed.
% If the numSlices argument is <=0, no data are read; if it's negative, the file
% is left open and the handle fh is returned instead of m.  Subsequent calls to fread() will then load the
% data.  For example,
%   [h s]=ReadMRC(name,1,-1);
%   for i=1:s.nz
%       m1=fread(h,s.nx*s.ny,s.string);
%       m1=reshape(m1,s.nx,s.ny);
%        --do something with m1--
%   end;
%   fclose(h);
% If a third output argument is specified, the header (as a vector of
% int32) is returned.  For non-extended files it has 256 entries.
%
% Added interpretation of the extended header (length = c(2) which is important with imod and
% serialEM at Brandeis.  lw 15 Aug 07
% Added error tests.  fs 2 Sep 07
% Added the ability to read selected slices, and read big-endian files
% fs 17 Sep 09
% Changed the a array so that the returned s values are doubles.  5 Nov 09
% Added data mode 6 (unsigned int16) fs 20 Jan 10
% Changed s fields to doubles fs 20 Feb 10
% Added fields s.mx, s.my and s.mz (unit cell dimensions), s.org (origin
% offset) along with s.pixA the pixel (voxel) size, which is computed as
% s.pixA=s.rez/s.mx.  The minimum, maximum and average pixel values are
% returned as s.mi, s.ma and s.ma.fs 10 Jul 11
%
% Added support for packed bytes (we called it mode 32 but SerialEM calls
% it mode 101. Both are recognized. The result is returned as uint8s fs 17 Sep 13 
% Added s.err to give an error message when the file is incomplete.
% s.err=1 means the header was read, but not enough data;
% s.err=2 means the file didn't even contain a complete header.

% Default returned values
hdr=[];
extraHeader=[];
s=struct;
s.err=0;
map=[];

if nargin<1 || numel(filename)<1 % no name given: put up a file selector
    fprintf('Select an MRC file: ');
    [filename,path]=uigetfile('*.mrc;*.mrcs','Open an MRC file');
    if isnumeric(filename)
        disp(' no file selected.');
        return
    else
        disp([path filename]);
    end;
    cd(path);
end;

if nargin<2 || startSlice<1
    startSlice=1;
end;
if nargin<3
    numSlices=inf;
end;
if nargin<4
    test=0;
end;


% We first try for little-endian data
f=fopen(filename,'r','ieee-le');
if f<0
    error(['in ReadMRC the file could not be opened: ' filename])
end;

% Check the file size
fseek(f,0,'eof');  % go to the end of the file
nBytes=ftell(f);
fseek(f,0,'bof');  % rewind to start of the file
if nBytes<1024  % Not enough room for the header
    s.err=2;  % Serious error
    warning('File is too short for mrc header');
    return
end;

% Get the first 10 values, which are integers:
% nc nr ns mode ncstart nrstart nsstart nx ny nz
a=fread(f,10,'*int32');

% Check the nx value.
if abs(a(1))>1e5  % we must have the wrong endian data.  Try again.
    fclose(f);
    f=fopen(filename,'r','ieee-be');
    a=fread(f,10,'int32');  % convert to doubles
end;

if test
    a(1:10)
end;

s.mode=a(4);

% Get the next 12 (entries 11 to 23), which are floats.
% the first three are the cell dimensions.
% xlength ylength zlength alpha beta gamma mapc mapr maps
% amin amax amean.
[b,cnt]=fread(f,12,'float32');
if test
    b
end;

% b(4,5,6) angles
s.mi=b(10); % minimum value
s.ma=b(11); % maximum value
s.av=b(12);  % average value

s.rez=double(b(1)); % cell size x, in A.

% get the next 30, which brings us up to entry no. 52.
[c,cnt]=fread(f,30,'*int32');
if test, c(1:3), end;
% c(2) is the extended header in bytes.

% the next two are supposed to be character strings.
[d,cnt]=fread(f,8,'*uint8');
s.chars=char(d)';
if test
    d
end;

% Two more ints...
[e,cnt]=fread(f,2,'*int32');
if test, e, end;

% up to 10 strings....
ns=min(e(2),10);
str=char(zeros(10,80));
for i=1:10
    [g,cnt]=fread(f,80,'*uint8');
    str(i,:)=char(g)';
end;
% disp('header:'); disp(' ');
% disp(str(1:ns,:));
% disp(' ');
s.header=str(1:ns,:);

% Get ready to read the data.
s.nx=double(a(1));
s.ny=double(a(2));
s.nz=double(a(3));
s.mx=double(a(8));
s.my=double(a(9));
s.mz=double(a(10));
s.org=double(a(5:7));
s.pixA=s.rez/s.mx;

nx1=s.nx;
packedBytes=0;
switch s.mode
    case 0
        string='*uint8';
        pixbytes=1;
    case 1
        string='*int16';
        pixbytes=2;
    case 2
        string='*float32';
        pixbytes=4;
    case 6
        string='*uint16';
    case {32 101}
        packedBytes=1;
        string='*uint8';
        pixbytes=1;
        nx1=ceil(nx1/2);  % no. of bytes to read per row
    otherwise
        error(['ReadMRC: unknown data mode: ' num2str(s.mode)]);
end;

s.string=string;

% Make sure we are at the end of the header.
fseek(f,1024,'bof');

if nBytes<1024+c(2)  % not enough space for extended header
    s.err=2;
    warning('File is too short for the mrc extended header');
    return
end;

if(c(2)>0)
%     [extraHeader,cnt]=fread(f,c(2)/4,'*float32');
    [extraHeader,cnt]=fread(f,c(2)/2,'*int16');
    if test
%         disp(['Read extra header of ',num2str(c(2)/4),' floats!'])
        disp(['Read extra header of ',num2str(c(2)/2),' int16s'])
        extraHeader(1:16)
    end;
    %    disp((ex_header'));
end

if startSlice>1
    skipbytes=(startSlice-1)*nx1*s.ny*pixbytes;
    fseek(f,skipbytes,'cof');
end;

if numSlices<0
    map=f;
    return;
elseif numSlices==0
    fclose(f);
    map=0;
    return
end;

nz=double(max(0,min(s.nz-(startSlice-1),numSlices)));
% Must be double because ndata can require double precision
nDataToRead=nx1 * s.ny * nz;

if test
    string
    nDataToRead
end;

[map,dcnt]=fread(f,nDataToRead,string);

if dcnt ~= nDataToRead
    warning(sprintf('Not enough data in file: %d pixels, %d expected.',dcnt,nDataToRead));
    s.err=1;
end;
map(dcnt+1:nDataToRead)=0;  % fill up the rest of the data array

if packedBytes  % have to unpack bytes
    tempMap=rem(map',16);
    tempMap(2,:)=idivide(map',16);
    if 2*nx1>s.nx  % there's a padded element
%         diff= 2*nx1-s.nx
        map=reshape(tempMap,[nx1*2 s.ny nz]);
        map(nx1*2,:,:)=[]; % delete one element
    else
        map=reshape(tempMap,[s.nx s.ny nz]);
    end;
end;

if nargout>2 % We read an image of the header
    frewind(f);
    [hdr, cnt]=fread(f,256+c(2)/4,'*int32');
end;

fclose(f);

map=reshape(map,s.nx,s.ny,nz);

%
% f=fopen(filename,'r','ieee-be');
% fseek(f,1024,'bof');
%
% if(c(2)>0)
%     [extraHeader,cnt]=fread(f,c(2)/4,'*float32');
%     disp(['Read extra header of ',num2str(c(2)/4),' floats!'])
% %    disp((ex_header'));
% end
% fclose(f);