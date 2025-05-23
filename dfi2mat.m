function [output, debug] = dfi2mat(inFile)
%DFI2MAT Loads a dfi file into a MATLAB structure
%
%   [output, debug] = dfi2mat(inFile)
%
% dfi2mat loads the dfi specified in the string INFILE and attempts to
% convert it into a MATLAB structure.  OUTPUT may contain some combination
% of the following fields:
%   -- imageList    a list of all image planes found.
%   -- (image)      an array containing the image named (image), where
%                   (image) is an element of imageList.
%   -- colormap     a colormap for the dfi for use with the MATLAB function
%                   colormap.
%   -- vectors      a cell array containing the components of a set of
%                   vectors. Vectors are named as in the original dfi. If
%                   rescaling has been used then vectors{k}.x and
%                   vectors{k}.y contain positions for use with the MATLAB
%                   function quiver.
% Other fields may also be present but are less likely to be of use.
%
% DEBUG contains a complete record of everything found in the dfi file
% without any significant processing. If DFI2MAT gives any messages about
% duplicate objects then DEBUG will still contain the data. DEBUG contains
% the following fields
%   -- dfi2mat      structure containing fields
%       -- version  version of dfi2mat used
%       -- date     date of dfi2mat release
%   -- idFormat     string used in identification of dfi files.
%   -- Version      integer used in identification of dfi files.
%   -- Tag{k}       cell array of data tags. Contains fields
%       -- DataType integer containing details of the data type in object k
%       -- numBytes size of Object{k} in bytes. Don't trust.
%   -- Object{:}    cell array of object details. See Digiflow manual for
%                   meanings of the fields.
%
% DFI2MAT has not been tested on many types of dfi and is likely to fail
% for many of them. Not all features have been fully implemented and seeing
% warnings on running is normal. If/when you find a bug please send me an
% email at jsc63@cam.ac.uk.

if verLessThan('matlab','8.3')
    warning(['MATLAB version >= 8.3 (R2014a) is required for certain ' ...
        'functions. Consider upgrading your installation.']);
end

narginchk(1,1); % Need exactly one input

debug.dfi2mat.version = '1.6';
debug.dfi2mat.date = '2015/06/04';

if ~exist(inFile,'file');
    error([inFile ' does not exist']);
end

f = fopen(inFile,'r','l'); % Read as little-endian

if f == -1
    error(['Could not open ' inFile]);
end

% Read header
debug.idFormat = fread(f,32,'*char')';
debug.Version = fread(f,1,'*int32');

% Check that we've actually got a dfi file.
if ~strcmp(debug.idFormat,'Tagged floating point image file') ...
        || debug.Version~=0
    fclose(f);
    error([inFile ' does not appear to be a valid dfi file']);
end

% Used to warn if multiple tags of the same class have been found. This
% probably shouldn't happen but it pays to be safe.
imageFound = false;
range = false;
rescale = false;
colourScheme = false;
description = false;
userComments = false;
creatingProcess = false;
creatorDetails = false;
imageTime = false;
imageCoordinates = false;
imagePlaneDetails = false;

% Prep outputs in case things go wrong
output = [];
debug = [];

% MATLAB doesn't have do-while loops so we fake it
i = 1;
debug.Tag{1}.DataType = fread(f,1,'*int32');
debug.Tag{1}.nBytes = fread(f,1,'*int32'); % Don't trust this number


while(~feof(f)) % Until we run out of bytes
    
    % Object is constructed in a temporary structure o for clarity
    
    dataType = debug.Tag{i}.DataType; % For clarity
    
    switch dataType
        
        % Single plane image
        case {hex2dec('1001'), hex2dec('1004'), hex2dec('1008')}
            
            o.nx = fread(f,1,'*int32');
            o.ny = fread(f,1,'*int32');
            
            switch dataType
                case hex2dec('1001')
                    o.Type = '8 bit image';
                    o.r = fread(f,nx*o.ny,'int8=>double');
                    
                case hex2dec('1004')
                    o.Type = '32 bit image';
                    o.r = fread(f,o.nx*o.ny,'*single');
                    
                case hex2dec('1008')
                    o.Type = '64 bit image';
                    o.r = fread(f,o.nx*o.ny,'*double');
            end
            
            % We then reshape to give a standard MATLAB image - first two
            % indices are location, third selects colour. Finally we rotate
            % the matrix as DigiFlow uses a different index convention to
            % MATLAB.
            o.Image = rot90(reshape(o.r,[o.nx, o.ny]));
            
            if imageFound
                warning('Multiple image objects found. Ignoring.');
            else
                output.image = o.Image;
                output.imageType = o.Type;
                imageFound=true;
            end
            
            % Multi-plane images
        case {hex2dec('11001'), hex2dec('11004'), hex2dec('11008')}
            
            o.nx = fread(f,1,'*int32');
            o.ny = fread(f,1,'*int32');
            o.nz = fread(f,1,'*int32');
            
            switch dataType
                case hex2dec('11001')
                    o.Type = '8 bit multi-plane image';
                    o.r = fread(f,o.nx*o.ny*o.nz,'int8=>double');
                    
                case hex2dec('11004')
                    o.Type = '32 bit multi-plane image';
                    o.r = fread(f,o.nx*o.ny*o.nz,'*single');
                    
                case hex2dec('11008')
                    o.Type = '64 bit multi-plane image';
                    o.r = fread(f,o.nx*o.ny*o.nz,'*double');
            end
            
            % Convert to MATLAB image. See single plane image for details.
            tmp=reshape(o.r,[o.nx, o.ny, o.nz]);
            tmp2=zeros(o.ny,o.nx,o.nz);
            for k=1:size(tmp,3)
               tmp2(:,:,k) = rot90(squeeze(tmp(:,:,k)));
            end
            o.Image = tmp2;
            clear tmp tmp2
            
            if imageFound
                warning('Multiple image objects found. Ignoring.');
            else
                output.image = o.Image;
                output.imageType = o.Type;
                imageFound=true;
            end
            
            % Compressed images
        case {hex2dec('12001'), hex2dec('12004'), hex2dec('12008')}
            
            o.nx = fread(f,1,'*int32');
            o.ny = fread(f,1,'*int32');
            o.nz = fread(f,1,'*int32');
            
            o.szCompressed = fread(f,1,'*int32');
            o.c = fread(f,o.szCompressed,'*uint8');
            
            switch dataType
                case hex2dec('12001')
                    o.Type = 'Compressed 8bit image';
                    o.d = zlibdecode(o.c,'int8'); % Decompress and cast to
                    % correct type
                case hex2dec('12004')
                    o.Type = 'Compressed 32 bit image';
                    o.d = zlibdecode(o.c,'single');
                    
                case hex2dec('12008')
                    o.Type = 'Compressed 64 bit image';
                    o.d = zlibdecode(o.c,'double');
            end
            
            % Convert to MATLAB image. See single plane image for details.
            o.Image = rot90(reshape(o.d,[o.nx o.ny o.nz]));
            
            if imageFound
                warning('Multiple image objects found. Ignoring.');
            else
                output.image = o.Image;
                output.imageType = o.Type;
                imageFound=true;
            end
            
            % Image ranges
        case {hex2dec('1014'), hex2dec('1018')}
            switch dataType
                case hex2dec('1014')
                    o.Type = '32 bit range';
                    o.rBlack = fread(f,1,'*single');
                    o.rWhite = fread(f,1,'*single');
                case hex2dec('1018')
                    o.Type = '64 bit range';
                    o.rBlack = fread(f,1,'*double');
                    o.rWhite = fread(f,1,'*double');
            end
            
            if range
                warning('Multiple range objects found. Ignoring.')
            else
                output.imageRange = [o.rBlack o.rWhite];
                range = true;
            end
            
            % Rescalings
        case hex2dec('1100')
            o.Type = 'Rescale image';
            
            o.nxWant = fread(f,1,'*int32');
            o.nyWant = fread(f,1,'*int32');
            o.method = fread(f,1,'*int32');
            
            if rescale
                warning('Multiple rescaling objects found. Ignoring.')
            else
                output.want = double([o.nxWant o.nyWant]);
                output.method = o.method;
                output.useRectangle = false;
                rescale = true;
            end
            
        case hex2dec('1101')
            o.Type = 'Rescale image rectangle';
            
            o.nxWant = fread(f,1,'*int32');
            o.nyWant = fread(f,1,'*int32');
            o.method = fread(f,1,'*int32');
            o.useRectangle = fread(f,1,'*int32');
            o.Rectangle.xMin = fread(f,1,'*int32');
            o.Rectangle.yMin = fread(f,1,'*int32');
            o.Rectangle.xMax = fread(f,1,'*int32');
            o.Rectangle.yMax = fread(f,1,'*int32');
            
            if rescale
                warning('Multiple rescaling objects found. Ignoring.')
            else
                output.want = double([o.nxWant o.nyWant]);
                
                output.method = o.method;
                output.useRectangle = o.useRectangle;
                if(output.useRectangle)
                    % Use MATLAB style rectangle variables
                    % [x0 y0 width height]
                    output.rectangle = double([ ...
                        o.Rectangle.xMin, ...
                        o.Rectangle.yMin, ...
                        o.Rectangle.xMax - o.Rectangle.xMin, ...
                        o.Rectangle.yMax - o.Rectangle.yMin ]);
                end
                rescale = true;
                
            end
            
            % Colour schemes
        case hex2dec('2000')
            o.Type = 'Colour scheme';
            
            o.red = fread(f,256,'uint8');
            o.green = fread(f,256,'uint8');
            o.blue = fread(f,256,'uint8');
            
            if colourScheme
                warning('Multiple colour schemes detected. Ignoring.');
            else
                output.colormap = double([o.red o.green o.blue]) ./ 255;
                output.colormapName = 'DFI specified';
                colourScheme = true;
            end
            
        case hex2dec('2001')
            o.Type = 'Colour scheme name';
            
            o.name = fread(f,64,'*char')';
            
            if colourScheme
                warning('Multiple colour schemes detected. Ignoring.');
            else
                output.colormap = digiflowColourScheme(o.name);
                if ~isempty(output.colormap)
                    output.colormapName = o.name;
                    colourScheme = true;
                else
                    output = rmfield(output,'colormap');
                end
            end
            
        case hex2dec('2002')
            o.Type = 'Colour scheme name variable';
            
            o.iLen = fread(f,1,'*int32');
            o.name = fread(f,o.iLen,'*char')';
            
            if colourScheme
                warning('Multiple colour schemes detected. Ignoring.');
            else
                output.colormap = digiflowColourScheme(char(o.name));
                if ~isempty(output.colormap)
                    output.colormapName = o.name;
                    colourScheme = true;
                else
                    output = rmfield(output,'colormap');
                end
            end
            
            % Descriptive
        case hex2dec('3000')
            o.Type = 'Description';
            
            o.Descr = fread(f,512,'*char')';
            
            if description
                warning('Multiple descriptions detected. Ignoring.');
            else
                output.description = o.Descr;
                description = true;
            end
            
        case hex2dec('3001')
            o.Type = 'User comments';
            
            o.nBytes = fread(f,1,'*int32');
            o.Descr = fread(f,o.nBytes,'*char')';
            
            if userComments
                warning('Multiple user comments detected. Ignoring.');
            else
                output.userComments = o.Descr;
                userComments = true;
            end
            
        case hex2dec('3002')
            o.Type = 'Creating process';
            
            o.nBytes = fread(f,1,'*int32');
            o.Descr = fread(f,o.nBytes,'*char')';
            
            if creatingProcess
                warning('Multiple creating processes found. Ignoring.');
            else
                output.creatingProcess = o.Descr;
                creatingProcess = true;
            end
            
        case hex2dec('3003')
            o.Type = 'Creator details';
            
            o.Digiflow = fread(f, 32, '*char')';
            o.buildDate = fread(f,16,'*char')';
            o.licenceType = fread(f,16,'*char')';
            o.nameUser = fread(f,32,'*char')';
            o.nameComputer = fread(f,32,'*char')';
            o.nameDomain = fread(f,32,'*char')';
            o.guidUser = fread(f,32,'uint8=>char')';
            
            o.macAddress = fread(f,[12,4],'*char')';
            o.ipAddress = fread(f,[16,4],'*char')';
            
            if creatorDetails
                warning(['Multiple sets of creator details found. ' ...
                    'Ignoring.']);
            else
                output.creatorDetails.digiflow = o.Digiflow;
                output.creatorDetails.buildDate = o.buildDate;
                output.creatorDetails.licenceType = o.licenceType;
                output.creatorDetails.userName = o.nameUser;
                output.creatorDetails.computerName = o.nameComputer;
                output.creatorDetails.domainName = o.nameDomain;
                output.creatorDetails.userGUID = o.guidUser;
                output.creatorDetails.MACaddress = o.macAddress;
                output.creatorDetails.IPaddress = o.ipAddress;
                
                creatorDetails = true;
            end
            
        case hex2dec('3018')
            o.Type = 'Image time';
            
            o.iFrame = fread(f,1,'*int32');
            o.Reserved = fread(f,1,'*int32');
            o.Time = fread(f,1,'*double');
            o.tStep = fread(f,1,'*double');
            o.tFirst = fread(f,1,'*double');
            
            if imageTime
                warning('Multiple image time objects found. Ignoring.');
            else
                output.n = o.iFrame;
                output.t = o.Time;
                output.dt = o.tStep;
                output.t0 = o.tFirst;
                
                imageTime = true;
            end
            
            % Image details
        case hex2dec('4008')
            o.Type = 'Image coordinates';
            
            o.Kind = fread(f,1,'*int32');
            o.xWorldPerPixel = fread(f,1,'*double');
            o.yWorldPerPixel = fread(f,1,'*double');
            o.xOriginWorld = fread(f,1,'*double');
            o.yOriginWorld = fread(f,1,'*double');
            o.xUnits = fread(f,16,'*char')';
            o.yUnits = fread(f,16,'*char')';
            o.OriginalName = fread(f,64,'*char')';
            
            if imageCoordinates
                warning('Multiple coordinate systems detected. Ignoring.');
            else
                output.coordinates.kind = o.Kind;
                output.coordinates.xWorldPerPixel = o.xWorldPerPixel;
                output.coordinates.yWorldPerPixel = o.yWorldPerPixel;
                output.coordinates.xOriginWorld = o.xOriginWorld;
                output.coordinates.yOriginWorld = o.yOriginWorld;
                output.coordinates.xUnits = o.xUnits;
                output.coordinates.yUnits = o.yUnits;
                output.coordinates.name = o.OriginalName;
                
                imageCoordinates = true;
            end
            
        case hex2dec('4108')
            o.Type = 'Image plane details';
            
            o.nPlanes = fread(f,1,'*int32');
            
            % Preallocate
            o.Contains=zeros(o.nPlanes,1,'int32');
            o.Descr = cell(o.nPlanes,1);
            o.ParamA = zeros(o.nPlanes,1,'double');
            o.ParamB = zeros(o.nPlanes,1,'double');
            o.ParamC = zeros(o.nPlanes,1,'double');
            o.ParamD = zeros(o.nPlanes,1,'double');
            o.Unknown = cell(o.nPlanes,1); % Not specified in spec.
                                           % Colour scheme?
            
            for j=1:o.nPlanes
                o.Contains(j) = fread(f,1,'*int32');
                o.Descr(j) = {fread(f,32,'*char')'};
                o.ParamA(j) = fread(f,1,'*double');
                o.ParamB(j) = fread(f,1,'*double');
                o.ParamC(j) = fread(f,1,'*double');
                o.ParamD(j) = fread(f,1,'*double');
                o.Unknown(j) = {fread(f,32,'*char')'};
            end
            
            if imagePlaneDetails
                warning(['Multiple image plane details objects found.' ...
                    ' Ignoring.']);
            else
                output.plane=cell(o.nPlanes,1);
                
                for j=1:o.nPlanes
                    output.plane{j}.type = o.Contains(j);
                    output.plane{j}.name = char(strtrim(o.Descr(j)));
                    output.plane{j}.a = o.ParamA(j);
                    output.plane{j}.b = o.ParamB(j);
                    output.plane{j}.c = o.ParamC(j);
                    output.plane{j}.d = o.ParamD(j);
                    output.plane{j}.unknown = char(o.Unknown(j));
                end
                
                imagePlaneDetails = true;
            end
            
        otherwise
            if debug.Tag{i}.DataType < 0
                warning(['Negative data type: -#' ...
                    dec2hex(-debug.Tag{i}.DataType) '.']);
                warning('Later objects and tags are probably corrupt.');
            else
                warning(['Unknown data type: #' ...
                    dec2hex(debug.Tag{i}.DataType) '.']);
            end
            
            o.Raw = fread(f,debug.Tag{i}.nBytes,'*char');
    end
    
    % Store o and clear it for the next iteration
    debug.Object{i} = o;
    clear('o');
    
    i = i+1;
    debug.Tag{i}.DataType = fread(f,1,'*int32');
    debug.Tag{i}.nBytes = fread(f,1,'*int32');
end

% If we hit eof then the last tag is junk
debug.Tag(i) = [];

% Close the file
fclose(f);


% Now for some processing

if imageFound && ~isfield(output,'plane') % Single plane image
    output.imageList = {'image'};
    
elseif isfield(output,'plane') % Multi-plane images need some assembly work
    raw = output.image;
    output = rmfield(output,'image');
    
    numImages = 0;
    numGrey = 0;
    numRed = 0;
    numBlue = 0;
    numGreen = 0;
    
    numX = 0;
    numY = 0;
    numZ = 0;
    
    numU = 0;
    numV = 0;
    numW = 0;
    
    for k=1:length(output.plane)
        switch output.plane{k}.type
            case hex2dec('000')
                warning(['Plane ' num2str(k) ' has no type.']);
                
            case {hex2dec('001') hex2dec('010')}
                % Type #001 is a greyscale image
                
                % Type #010 is unspecified but is probably a scalar field
                % so treating as greyscale image
                
                if isfield(output,output.plane{k}.name)
                    warning([ ...
                        'Multiple images with the same name detected. ' ...
                        'Ignoring.']);
                else
                    numImages = numImages + 1;
                    numGrey = numGrey + 1;
                    validatedName = ...
                        makeValidFieldName(output.plane{k}.name,output);
                    output.imageList{numImages} = validatedName;
                    output.(validatedName) = raw(:,:,k);
                    
                    % Not in spec but probably correct
                    output.colormap = ...
                        digiflowColourScheme(output.plane{k}.unknown);
                end
                
            case hex2dec('002')
                numRed = numRed + 1;
                red{numRed}.im = raw(:,:,k);
                red{numRed}.name = (output.plane{k}.name);
                
            case hex2dec('003')
                numGreen = numGreen + 1;
                green{numGreen}.im = raw(:,:,k);
                green{numGreen}.name = (output.plane{k}.name);
                
            case hex2dec('004')
                numBlue = numBlue + 1;
                blue{numBlue}.im = raw(:,:,k);
                blue{numBlue}.name = (output.plane{k}.name);
                
            case hex2dec('101')
                numX = numX + 1;
                coords{numX}.x  = raw(:,:,k);
                coords{numX}.xname = (output.plane{k}.name);
                
            case hex2dec('102')
                numY = numY + 1;
                coords{numY}.y  = raw(:,:,k);
                coords{numY}.yname = (output.plane{k}.name);
                
            case hex2dec('103')
                numZ = numZ + 1;
                coords{numZ}.z  = raw(:,:,k);
                coords{numZ}.zname = (output.plane{k}.name);
                
            case hex2dec('201')
                numU = numU + 1;
                vectors{numU}.u = raw(:,:,k);
                vectors{numU}.uname = (output.plane{k}.name);
                
            case hex2dec('202')
                numV = numV + 1;
                % Take negative due to differences in index convention
                vectors{numV}.v = -raw(:,:,k);
                vectors{numV}.vname = (output.plane{k}.name);
                
            case hex2dec('203')
                numW = numW + 1;
                vectors{numW}.w = -raw(:,:,k);
                vectors{numW}.wname = (output.plane{k}.name);
                
            otherwise
                warning(['Unknown plane type #' ...
                    dec2hex(output.plane{k}.type) ...
                    '. Ignoring.']);
        end
    end
    
    % Reassemble images
    if (numRed ~= numBlue || numBlue ~= numGreen) % I.e. missing colours
        warning(['Incomplete coloured images detected.'...
            ' Ignoring incomplete images.']);
        
    end
    
    numColoured = min([numRed numBlue numGreen]);
    
    for k = 1:numColoured
        numImages = numImages + 1;
        validatedName = makeValidFieldName('image',output);
        output.imageList{numImages} = validatedName;
        output.(validatedName) = cat(3,red{k}.im,green{k}.im,blue{k}.im);
        output.componentNames{k}.r = red{k}.name;
        output.componentNames{k}.g = green{k}.name;
        output.componentNames{k}.b = blue{k}.name;
    end
    
    % Only allow all 3D or all 2D for simplicity
    if (numX ~= numY ||... %Unequal X and Y
            (numZ ~= 0 && numX ~= numZ)) % Or have 3D
                                         % and unequal X, Y and Z
        warning(['Incomplete coordinates detected.'...
            ' Ignoring incomplete coordinates.']);
        
    else
        if numZ > 0
            numCoords = min([numX numY numZ]);
        else
            numCoords = min([numX numY]);
        end
                
        % This orders the output of fieldnames(output.coordinates{k}) such
        % that it gives the components as {x, y, z}
        
        for k=1:numCoords
            names = ... 
                makeValidFieldName({coords{k}.xname, coords{k}.yname});
            output.coordinates{k}.(names{1}) = coords{k}.x;
            output.coordinates{k}.(names{2}) = coords{k}.y;
            
            if numZ > 0
                validatedZName = ...
                    makeValidFieldName(coords{k}.zname,output.coordinates{k});
                output.coordinates{k}.(validatedZName) = coords{k}.z;
            end
            
        end
    end
    
    % Reassemble vectors
    % Only allow all 3D or 2D
    if (numU ~= numV ||... %Unequal U and V
            (numW ~= 0 && numU ~= numW)) % Or have 3D and unequal U, V and W
        warning(['Incomplete vector fields detected.'...
            ' Ignoring incomplete vector fields.']);
        
    else
        if numW > 0
            numVectors = min([numU numV numW]);
        else
            numVectors = min([numU numV]);
        end
        
        % As with coordinates
        for k=1:numVectors
            names = ...
                makeValidFieldName({vectors{k}.uname,vectors{k}.vname});
            output.vectors{k}.(names{1}) = vectors{k}.u;
            output.vectors{k}.(names{2}) = vectors{k}.v;
            
            if numZ > 0                
                validatedWName = ...
                    makeValidFieldName(vectors{k}.wname,output.vectors{k});
                output.vectors{k}.(validatedWName) = vectors{k}.w;
            end
        end
    end
end

if rescale
    
    if output.useRectangle
        rectangle = output.rectangle;
    else
        rectangle = [1, 1, output.want(1), output.want(2)];
    end
    
    width = rectangle(3) - rectangle(1);
    height = rectangle(4) - rectangle(2);
    
    % Find closest MATLAB interpolation method.
    switch output.method
        case 0 % Constant
            interpType = 'nearest';
        case 1 % Bilinear
            interpType = 'linear';
        case 2 % Bicubic
            interpType = 'cubic';
        case 3 % Natural spline
            %warning('Natural splines unsupported. Using MATLAB splines.');
            interpType = 'spline';
        case 4 % Cubic b-spline
            %warning('Cubic b-splines unsupported. Using MATLAB splines.');
            interpType = 'spline';
        case 5 % Quintic b-spline
            %warning('Quintic b-splines unsupported. Using MATLAB splines.');
            interpType = 'spline';
        otherwise
            warning('Unknown rescaling type. Using linear interpolation.')
            interpType = 'linear';
    end
    
    if isfield(output,'imageList')
        for j=1:length(output.imageList)
            % imageList should be valid already
            im = output.imageList{j};
            % MATLAB size function gives number of rows by number of columns
            % I.e. y resolution by x resolution
            [y, x, z] = size(output.(im));
            
            X = 1:x;
            Y = 1:y;
            
            XX = 1:((x-1)/width):x;
            YY = 1:((y-1)/height):y;
            
            old = output.(im);
            output.(im) = zeros(output.want(2),output.want(1),z);
            
            for k=1:z % Loop over components
                output.(im)((rectangle(2):rectangle(4))+1, ...
                (rectangle(1):rectangle(3))+1, k) ...
                    = interp2(X, Y, old(:,:,k), XX, YY', interpType);
            end
        end
    end
    
    if isfield(output,'vectors')
        for j=1:length(output.vectors)
            
            components = fieldnames(output.vectors{j});
                       
            [y, x, z] = size(output.vectors{j}.(components{1}));
            
            if z ~= 1
                warning(['Cannot deal with 3D vector rescaling - ' ...
                    'rescale manually.']);
                continue;
            end
            
            X = 1:x;
            Y = 1:y;
            
            XX = 1:((x-1)/width):x;
            YY = 1:((y-1)/height):y;        
            
            
            for k=1:length(components)
                temp = zeros(output.want(2), output.want(1));
                temp((rectangle(2):rectangle(4))+1, (rectangle(1):rectangle(3))+1) ...
                    = interp2(X, Y, output.vectors{j}.(components{k}), XX, YY', interpType);
                output.vectors{j}.(components{k}) = temp;
            end
            
        end
    end
end

% If we've only got one vector field then copy into main structure for ease
% of use

if isfield(output,'vectors') && (length(output.vectors) == 1)
    vectors = output.vectors{1};
    components = fieldnames(vectors);
    % Check for clashes in field names
    if ~any(isfield(output,components)) 
        for k=1:length(components)
            output.(components{k}) = vectors.(components{k});
        end
    end
end

end

function output = digiflowColourScheme(name)
%DIGIFLOWCOLOURSCHEME Implements standard Digiflow colour schemes
%
%   output = digiflowColourScheme(name)
%
% NAME must be a character array.  OUTPUT is a valid MATLAB colormap.

% Ignore trailing spaces and nulls that can appear
switch deblank(name)
    
    case 'random - random'
        red=[ 0.7804 0.6275 0.7216 0.3255 0.9804 0.8118 0.3569 0.5059 0.6118 0.8706 0.7608 0.2157 0.0000 0.4275 0.7451 0.4000 0.8902 0.2471 0.2902 0.2353 0.9451 0.8863 0.7059 0.9255 0.1922 0.4510 0.0157 0.0353 0.1608 0.6510 0.7059 0.4706 0.8118 0.3020 0.5804 0.4706 0.4275 0.4353 0.9608 0.7961 0.0706 0.5804 0.7412 0.3373 0.2275 0.5020 0.0275 0.6902 0.3804 0.7608 0.8000 0.9451 0.1608 0.5255 0.5961 0.2157 0.6510 0.0863 0.8353 0.1059 0.4353 0.6667 0.0314 0.5608 0.6000 0.7020 0.6510 0.2706 0.7804 0.5804 0.8157 0.1569 0.1451 0.9608 0.4000 0.4863 0.1804 0.0000 0.4706 0.1255 0.8000 0.3804 0.2824 0.0667 0.1255 0.8902 0.3412 0.1608 0.3451 0.7765 0.4510 0.8275 0.7059 0.1569 0.6706 0.3216 0.0118 0.6000 0.2824 0.6157 0.1765 0.8510 0.9255 0.3451 0.8902 0.5451 0.4353 0.1608 0.7608 0.1451 0.6706 0.3255 0.0000 0.4118 0.6902 0.5608 0.6510 0.1216 0.9569 0.8667 0.6353 0.6824 0.7804 0.2510 0.3216 0.8314 0.8118 0.5412 0.6902 0.3961 0.1804 0.3804 0.0902 0.2902 0.6353 0.0000 0.1608 0.4510 0.6706 0.0118 0.1961 0.2000 0.7059 0.9804 0.2000 0.6275 0.7569 0.7804 0.5804 0.1373 0.5059 0.6510 0.6314 0.3255 0.2510 0.6510 0.4157 0.8863 0.3608 0.2510 0.1765 0.9451 0.8353 0.1608 0.6510 0.8863 0.0706 0.6706 0.0353 0.5961 0.1961 0.9255 0.2275 0.7451 0.9451 0.5373 0.9255 0.2353 0.6667 0.8510 0.5255 0.2510 0.4353 0.6000 0.0157 0.7059 0.0353 0.4706 0.0157 0.6157 0.2706 0.4353 0.0353 0.4902 0.5608 0.4510 0.7412 0.1451 0.5255 0.0902 0.9961 0.6471 0.8353 0.9412 0.2353 0.9961 0.7608 0.8510 0.4353 0.4353 0.9765 0.5961 0.6902 0.7059 0.9020 0.6510 0.8902 0.9608 0.2824 0.2824 0.7451 0.1608 0.8471 0.4706 0.4000 0.1255 0.8157 0.1961 0.1255 0.1255 0.9608 0.0353 0.0510 0.0706 0.8353 0.3451 0.2706 0.2000 0.8902 0.3059 0.3608 0.4275 0.8706 0.2706 0.6706 0.1255 0.7922 0.0510 0.4510 0.1961 0.0275 0.8353 0.1608 0.1608 0.2314 0.3255];
        green=[ 0.1804 0.3569 0.8471 0.2902 0.2314 0.3961 0.6863 0.7059 0.6510 0.7765 0.3412 0.7569 0.8275 0.1020 0.2706 0.3961 0.9804 0.7804 0.0275 0.4902 0.9255 0.0510 0.5961 0.2510 0.0000 0.2353 0.4118 0.8353 0.0902 0.5216 0.7373 0.4471 0.7412 0.7804 0.5059 0.3255 0.9804 0.4353 0.9608 0.5059 0.3804 0.9608 0.1804 0.6863 0.4353 0.5451 0.5804 0.4157 0.9804 0.3608 0.6314 0.4863 0.9922 0.9255 0.1961 0.1451 0.3804 0.0706 0.2275 0.0510 0.6902 0.3412 0.9373 0.2706 0.3373 0.7608 0.1804 0.5255 0.4667 0.5059 0.7608 0.0510 0.2667 0.6510 0.7451 0.2353 0.0471 0.4510 0.9804 0.0510 0.4471 0.4667 0.5255 0.7569 0.1412 0.5059 0.8000 0.7451 0.5922 0.9412 0.8000 0.1255 0.9804 0.7412 0.8353 0.8275 0.3922 0.3804 0.5373 0.9451 0.2706 0.9922 0.4353 0.7255 0.0118 0.2353 0.3922 0.7059 0.3059 0.2275 0.4275 0.8000 0.5059 0.1059 0.2863 0.1255 0.3804 0.3020 0.1451 0.6000 0.5059 0.5412 0.1412 0.5608 0.5255 0.7255 0.1059 0.8118 0.1451 0.5059 0.3804 0.2706 0.2902 0.8000 0.9216 0.8314 0.8314 0.2157 0.0353 0.9373 0.3255 0.1255 0.9922 0.7373 0.1608 0.2902 0.0314 0.1216 0.5373 0.0157 0.8353 0.5804 0.8000 0.3569 0.8000 0.9804 0.4510 0.4000 0.6157 0.2706 0.6510 0.1804 0.1255 0.5608 0.0353 0.7608 0.9608 0.9373 0.5059 0.6353 0.2157 0.8667 0.4706 0.4157 0.4275 0.8157 0.1255 0.7255 0.4902 0.9255 0.4510 0.1373 0.8667 0.0706 0.9059 0.5961 0.6353 0.2157 0.7608 0.2706 0.3059 0.5804 0.0824 0.3961 0.0510 0.3255 0.8314 0.2510 0.6000 0.7804 0.7216 0.4706 0.7608 0.8510 0.6706 0.8706 0.9255 0.7216 0.8000 0.2863 0.6667 0.0353 0.4902 0.8471 0.9804 0.1804 0.6353 0.6667 0.8510 0.4510 0.6667 0.7373 0.2863 0.5608 0.3961 0.9608 0.3765 0.2314 0.5059 0.8510 0.4510 0.7804 0.4510 0.3059 0.3255 0.9059 0.0510 0.4353 0.5451 0.2824 0.9059 0.9961 0.2000 0.0000 0.8157 0.1961 0.0902 0.6902 0.1961 0.8510 0.5451 0.8706 0.3216 0.5961 0.8510 0.3255];
        blue=[ 0.4510 0.7059 0.7922 0.4510 0.4706 0.4863 0.0824 0.7059 0.7255 0.4157 0.4275 0.1059 0.6157 0.9922 0.4510 0.7804 0.8510 0.2000 0.7961 0.5020 0.6510 0.0000 0.1765 0.5804 0.7059 0.6706 0.1451 0.4706 0.0667 0.5608 0.1804 0.5216 0.0863 0.4706 0.2157 0.5059 0.9608 0.6000 0.9255 0.1216 0.7569 0.9608 0.4353 0.9608 0.9804 0.4824 0.2706 0.4706 0.6000 0.0902 0.8510 0.2706 0.1216 0.1412 0.1059 0.1608 0.2510 0.7059 0.2824 0.9255 0.0667 0.4000 0.4706 0.1373 0.3451 0.6000 0.6902 0.3608 0.8706 0.1255 0.9059 0.4000 0.7059 0.1059 0.2275 0.4471 0.9608 0.9020 0.1216 0.0824 0.2510 0.2314 0.0118 0.6353 0.6667 0.0667 0.7608 0.6824 0.8471 0.0000 0.2275 0.8314 0.0667 0.0667 0.9059 0.3255 0.3451 0.9608 0.8902 0.3804 0.2314 0.0902 0.0510 0.8000 0.8667 0.1059 0.2902 0.5059 0.0706 0.9608 0.8824 0.9216 0.8902 0.0824 0.0510 0.3961 0.1059 0.5451 0.6157 0.4706 0.3804 0.8157 0.3451 0.7059 0.6510 0.4706 0.6824 0.4157 0.3255 0.6471 0.5059 0.8706 0.1961 0.0314 0.9451 0.7020 0.0706 0.9059 0.5059 0.6667 0.7255 0.7922 0.0510 0.9765 0.8000 0.5608 0.2902 0.0353 0.4353 0.9412 0.0706 0.5922 0.8706 0.3451 0.3922 0.6353 0.4902 0.0157 0.7255 0.8000 0.7608 0.0157 0.2902 0.0902 0.6510 0.0000 0.9608 0.4902 0.7608 0.4510 0.0706 0.2353 0.3608 0.3804 0.3804 0.7922 0.6157 0.6353 0.8510 0.3412 0.1451 0.8275 0.2275 0.4706 0.2863 0.0706 0.3804 0.6863 0.5451 0.8902 0.3059 0.8314 0.0353 0.4706 0.4275 0.2157 0.5059 0.1569 0.9451 0.0902 0.1804 0.6863 0.3804 0.3569 0.1059 0.3451 0.2118 0.7216 0.9059 0.2157 0.2863 0.0275 0.7569 0.6510 0.7059 0.3961 0.2275 0.0275 0.9608 0.1373 0.2353 0.0353 0.9020 0.9373 0.0353 0.6353 0.1608 0.5608 0.5765 0.8353 0.2157 0.6824 0.0706 0.5373 0.7451 0.2157 0.3804 0.3922 0.3059 0.2510 0.4902 0.0353 0.2157 0.1765 0.4510 0.6706 0.4706 0.9608 0.5961 0.0706 0.8118 0.2000 0.9608 0.5922 0.0824 0.4706];
    case 'random - inverted bipolar'
        red=[ 0.0000 0.0000 0.0000 0.0157 0.0118 0.0275 0.0824 0.0510 0.0353 0.1059 0.0510 0.0863 0.1765 0.0000 0.0667 0.1216 0.3059 0.1569 0.0706 0.1059 0.1020 0.2510 0.0902 0.2157 0.1451 0.2000 0.3608 0.2510 0.1765 0.1608 0.0510 0.2902 0.1373 0.3804 0.4118 0.2706 0.4353 0.0157 0.0510 0.2000 0.2510 0.2157 0.2510 0.3765 0.0000 0.6353 0.1451 0.4510 0.1451 0.5059 0.7216 0.3412 0.8510 0.2157 0.5412 0.5255 0.0706 0.5804 0.5569 0.6667 0.5020 0.2706 0.4667 0.0157 0.6000 0.9255 0.7451 0.6118 1.0000 0.2353 0.3451 0.7255 0.5804 0.1059 0.7255 0.5922 0.7255 0.5255 0.8275 0.3059 0.7059 1.0000 0.1765 1.0000 0.7961 0.5412 0.8510 0.3216 0.9255 0.7451 0.0824 1.0000 1.0000 0.3804 1.0000 0.8706 0.1451 0.2902 1.0000 1.0000 0.8902 1.0000 1.0000 1.0000 1.0000 0.1961 0.4706 1.0000 0.3922 0.3059 1.0000 1.0000 1.0000 0.8353 0.3922 1.0000 0.4706 1.0000 1.0000 1.0000 1.0000 0.1451 1.0000 0.7451 0.4902 0.6157 1.0000 0.5608 0.0510 0.4471 0.9804 0.8000 0.5961 1.0000 0.8314 1.0000 0.3451 0.6275 0.9373 0.9961 0.4157 1.0000 0.7255 1.0000 0.9451 1.0000 0.8157 0.6510 0.1608 0.1059 1.0000 0.6157 1.0000 0.4510 0.3451 1.0000 0.8706 0.7059 0.7922 0.3059 0.4157 1.0000 1.0000 0.7765 1.0000 1.0000 0.2157 0.6667 0.3569 0.5451 0.4706 0.8471 0.6510 1.0000 0.7608 0.6000 0.4353 0.6510 0.9020 0.8000 0.8157 0.9961 0.2157 1.0000 0.9255 0.6706 0.1255 0.6706 0.4353 1.0000 0.2863 1.0000 0.6157 0.8471 0.4706 0.5255 0.5412 0.6353 0.7059 0.7020 0.3804 0.5255 0.7451 0.1608 0.1373 0.4000 0.3922 0.1961 0.8157 0.2000 0.4353 0.2863 0.4706 0.0824 0.3216 0.3451 0.1451 0.2706 0.2863 0.3804 0.1059 0.2902 0.0510 0.0353 0.2902 0.2353 0.2902 0.0000 0.1804 0.2000 0.0314 0.1804 0.2510 0.2275 0.2706 0.3059 0.0824 0.1804 0.1412 0.1608 0.1373 0.1059 0.2314 0.1059 0.1020 0.0118 0.0706 0.0510 0.0510 0.0157 0.0157 0.0157 0.0000 0.0000 0.0000 0.0000];
        green=[ 0.0000 0.0000 0.0157 0.0353 0.0314 0.0000 0.0157 0.0353 0.0471 0.0824 0.1020 0.1020 0.0000 0.0275 0.1804 0.0510 0.0471 0.0902 0.2118 0.1451 0.1059 0.0000 0.2510 0.2275 0.0000 0.0667 0.1608 0.1373 0.1255 0.2510 0.2000 0.0667 0.2902 0.1059 0.1255 0.4314 0.4000 0.3216 0.2706 0.4353 0.3412 0.5059 0.3804 0.3608 0.5020 0.1059 0.6706 0.4000 0.5059 0.3451 0.2706 0.3059 0.0471 0.2510 0.5059 0.0863 0.4706 0.3020 0.5059 0.3255 0.5255 0.7059 0.4000 0.7255 0.7412 0.1059 0.7765 0.5961 0.3451 0.1451 0.8157 0.9412 0.7059 0.7961 0.1059 0.1059 0.2471 0.5255 0.2824 0.3961 0.5804 0.4314 0.8510 0.6118 0.7059 0.7216 0.5059 0.9412 0.5804 0.4706 1.0000 0.6510 0.4902 0.7412 0.3373 0.5255 0.9255 1.0000 0.2157 0.0706 0.2667 0.5255 0.3216 0.7804 0.6118 1.0000 1.0000 0.1608 1.0000 0.6510 0.0902 0.8510 0.7804 1.0000 1.0000 0.7059 1.0000 0.0902 0.2667 0.8118 1.0000 1.0000 0.0157 1.0000 0.5608 1.0000 1.0000 0.4510 1.0000 1.0000 1.0000 1.0000 0.6353 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.0471 1.0000 0.9804 0.4471 0.2314 0.0510 0.9804 0.3255 0.8706 1.0000 1.0000 0.2157 1.0000 0.5961 1.0000 1.0000 0.5059 0.5373 0.2510 0.7569 1.0000 1.0000 0.3961 0.1451 1.0000 0.1216 0.2902 0.6314 1.0000 1.0000 0.2000 0.6510 0.8706 0.7059 0.3451 0.2510 0.1451 1.0000 0.5608 0.7765 0.6510 0.9059 0.6157 0.6510 0.2902 0.7216 0.7608 1.0000 0.5412 0.7059 0.2824 0.8000 0.1020 0.4000 0.1059 0.4314 0.3059 0.3255 0.3255 0.4824 0.3608 0.4902 0.7255 0.4667 0.7020 0.2902 0.5255 0.3804 0.5255 0.2510 0.4706 0.1569 0.1922 0.3059 0.0706 0.0902 0.1804 0.5059 0.3451 0.3373 0.1412 0.3059 0.4510 0.4157 0.0471 0.2157 0.2000 0.0706 0.4902 0.1451 0.2275 0.2902 0.0000 0.0510 0.1255 0.0275 0.1059 0.1804 0.1451 0.1608 0.0824 0.1569 0.1255 0.0157 0.0902 0.0667 0.0353 0.0667 0.0510 0.0902 0.0510 0.0510 0.0314 0.0157 0.0157 0.0000 0.0000];
        blue=[ 0.0000 0.0000 0.0157 0.0000 0.0353 0.0824 0.0353 0.0667 0.0902 0.0118 0.0667 0.0510 0.0902 0.2667 0.0706 0.1608 0.0157 0.1451 0.1255 0.1765 0.2471 0.2314 0.1451 0.0863 0.4000 0.3059 0.0824 0.2314 0.3451 0.2510 0.4353 0.3569 0.3216 0.2706 0.2510 0.1059 0.0000 0.5255 0.5569 0.2706 0.3451 0.2314 0.3255 0.2510 0.5255 0.2902 0.2510 0.2353 0.4510 0.2863 0.1608 0.5451 0.3059 0.7608 0.2000 0.6706 0.7608 0.4353 0.2902 0.3804 0.3765 0.4353 0.5804 0.7255 0.1451 0.4902 0.0157 0.3608 0.0902 1.0000 0.4706 0.0000 0.3922 0.8000 0.8902 1.0000 0.8000 0.7451 0.7216 1.0000 0.5804 0.4118 0.8902 0.0000 0.4510 0.7255 0.6471 0.7765 0.5451 0.8667 0.3412 0.0275 0.3451 1.0000 0.3804 0.8157 1.0000 0.8353 0.9608 1.0000 1.0000 0.7608 0.7804 0.5922 0.7608 1.0000 0.3255 1.0000 0.5922 1.0000 1.0000 0.4863 0.1059 0.5804 1.0000 0.9608 1.0000 1.0000 1.0000 0.8157 0.7059 0.5961 1.0000 1.0000 1.0000 0.7961 0.6000 1.0000 0.7059 0.8706 0.6118 0.9451 1.0000 0.7059 0.9804 0.2314 1.0000 0.6510 0.5020 1.0000 1.0000 0.4000 1.0000 1.0000 1.0000 0.5451 1.0000 0.9922 0.5059 1.0000 1.0000 0.7804 0.6510 0.8000 0.5059 0.7059 0.9059 1.0000 0.7255 0.4510 0.6706 0.4863 0.9765 0.1804 0.2000 0.1765 1.0000 0.2353 0.4824 1.0000 0.8510 0.2353 0.5765 0.2314 0.8706 1.0000 0.1255 0.6000 0.1255 0.3059 0.0275 0.1255 0.8353 0.1569 0.0118 0.1961 0.3059 0.3804 0.4275 0.0706 0.4353 0.0353 0.4471 0.4902 0.5216 0.5608 0.5059 0.3922 0.1373 0.2471 0.4157 0.0000 0.0275 0.3451 0.7608 0.2471 0.3608 0.4000 0.0157 0.3961 0.4510 0.5451 0.2157 0.8157 0.5451 0.4000 0.2510 0.2510 0.2353 0.3059 0.4000 0.0471 0.2863 0.6510 0.2000 0.2471 0.3059 0.1451 0.2902 0.1608 0.2510 0.3765 0.2157 0.1451 0.1922 0.0353 0.1608 0.0863 0.0902 0.1059 0.0510 0.0824 0.0510 0.0667 0.0824 0.1804 0.0667 0.0706 0.0118 0.0510 0.0353 0.0353 0.0353 0.0157 0.0000 0.0000];
    case 'random - bipolar'
        red=[ 1.0000 1.0000 0.3961 1.0000 0.7608 1.0000 1.0000 0.8902 0.2863 0.8510 1.0000 0.2902 0.9961 0.3255 0.7059 0.6706 1.0000 0.2902 1.0000 0.2706 1.0000 0.7451 0.7608 1.0000 1.0000 0.2000 0.0118 0.5608 0.1608 0.6667 0.9255 0.3059 0.7216 0.9059 1.0000 0.8510 1.0000 0.2667 0.8510 0.6706 0.0706 0.3412 1.0000 0.5255 0.6510 0.5255 0.3255 0.6902 0.6510 0.6157 0.4353 0.6275 0.2157 0.2510 0.2000 0.0706 0.8353 0.9569 0.8667 0.1608 0.1608 1.0000 0.2902 0.8824 0.4510 0.6510 0.2667 0.1804 0.3451 0.4157 0.4353 0.1255 0.5451 0.9059 0.5059 0.4275 0.6275 0.3451 0.2118 0.4510 0.4706 0.5059 0.1059 0.2000 0.5451 0.3608 0.1059 0.0510 0.3765 0.2353 0.0706 0.2510 0.3804 0.2353 0.3804 0.3961 0.3059 0.3804 0.2353 0.3451 0.0275 0.3922 0.2510 0.1961 0.1020 0.1922 0.0510 0.1608 0.1608 0.0118 0.1373 0.0706 0.0000 0.1765 0.0706 0.0510 0.1059 0.1020 0.1059 0.0863 0.0314 0.0353 0.0275 0.0353 0.0353 0.0314 0.0000 0.0000 0.0000 0.0157 0.0314 0.0157 0.0000 0.0314 0.0353 0.0902 0.1059 0.0902 0.1255 0.0000 0.0353 0.0510 0.0863 0.0510 0.0667 0.0510 0.1608 0.3216 0.1608 0.1804 0.1059 0.1804 0.2824 0.1961 0.2157 0.1608 0.0510 0.1804 0.1451 0.4000 0.1922 0.4863 0.3451 0.0471 0.4157 0.2902 0.1608 0.3451 0.3608 0.4353 0.4863 0.6471 0.5608 0.1765 0.7608 0.4510 0.3922 0.0471 0.1569 0.6706 0.6000 0.4902 0.2157 0.6510 0.5059 0.2863 0.3451 0.1922 0.0157 0.6510 0.4824 0.7451 0.9765 0.7059 1.0000 0.2353 0.5451 0.7608 0.8157 0.0314 0.4353 0.7804 0.4902 0.3216 0.7255 0.5804 0.2353 0.6510 0.6353 0.1608 0.4706 0.4353 0.3451 0.3373 0.7765 0.4157 0.4902 0.8863 0.9451 1.0000 0.3451 0.8000 0.4667 0.9255 1.0000 0.4510 1.0000 1.0000 0.8157 0.1608 1.0000 0.5608 0.7608 0.4510 0.9255 0.5922 0.2510 0.9373 0.8353 0.9608 1.0000 1.0000 1.0000 0.5059 0.9059 1.0000 1.0000 0.6275 0.7804 1.0000 0.6706 1.0000 0.3412 1.0000 0.3059 1.0000];
        green=[ 1.0000 0.6157 1.0000 0.1608 1.0000 1.0000 0.6275 0.9608 1.0000 0.8157 0.6706 0.9922 0.4118 1.0000 0.9373 0.6510 1.0000 1.0000 0.0157 1.0000 0.4157 1.0000 0.7216 0.3608 0.0000 1.0000 1.0000 0.9961 1.0000 1.0000 0.6000 0.2902 0.6902 1.0000 0.4902 0.4510 0.4824 1.0000 0.4157 1.0000 0.3569 0.5569 0.8157 1.0000 0.2353 0.8902 0.8353 1.0000 0.4510 0.4824 0.3451 0.3216 1.0000 0.9412 0.7608 0.8000 0.1059 0.0667 0.6706 0.0824 0.6353 0.1451 0.5451 0.4157 0.3765 0.1804 0.8902 0.3451 0.6510 0.6510 0.4118 0.8510 0.0667 0.3059 0.2706 0.3804 0.0118 0.1608 0.8902 0.2275 0.2471 0.0902 0.5765 0.5255 0.0353 0.2000 0.2706 0.2706 0.2510 0.3961 0.4471 0.3451 0.2863 0.3451 0.2667 0.1412 0.0510 0.0902 0.2863 0.2510 0.3961 0.2000 0.1922 0.1451 0.2314 0.1608 0.1608 0.2353 0.2706 0.1804 0.1059 0.1451 0.0510 0.1451 0.2353 0.2000 0.1020 0.0118 0.0000 0.0510 0.1216 0.0510 0.0000 0.0275 0.0275 0.0157 0.0000 0.0000 0.0000 0.0000 0.0000 0.0353 0.0667 0.0902 0.0157 0.0353 0.0706 0.0706 0.0706 0.2000 0.0706 0.0353 0.1451 0.1451 0.0510 0.1059 0.1216 0.0118 0.0510 0.0000 0.1804 0.0667 0.1922 0.2000 0.2706 0.3059 0.2157 0.2902 0.2510 0.0863 0.4353 0.0471 0.4000 0.4510 0.0000 0.2510 0.2902 0.2510 0.4510 0.4510 0.2706 0.2902 0.3765 0.5608 0.1608 0.5059 0.4706 0.5608 0.3451 0.4510 0.5804 0.0510 0.9451 0.0314 0.7373 1.0000 0.2157 0.5059 0.8157 0.0667 0.4902 0.3020 0.5059 0.6863 0.0275 1.0000 0.2510 0.8000 0.1961 0.9373 0.7451 0.7255 0.4471 0.9804 0.9961 0.8353 0.9059 0.9608 0.7412 0.8863 0.8157 0.7608 0.7451 0.4000 0.8000 1.0000 1.0000 1.0000 0.7608 0.3451 0.7765 1.0000 0.7059 0.2118 1.0000 0.7804 0.6353 0.1059 0.3765 1.0000 0.6157 0.9608 0.9059 0.9608 0.1608 1.0000 0.9922 0.6157 0.5059 0.8667 1.0000 1.0000 0.1020 1.0000 0.9608 1.0000 0.4706 1.0000 0.6353 1.0000 1.0000 0.5059 0.6353 0.4000 1.0000 1.0000];
        blue=[ 0.3059 1.0000 0.9804 1.0000 0.3765 0.3804 0.0706 0.9608 1.0000 1.0000 0.7059 1.0000 1.0000 0.7569 1.0000 1.0000 0.3255 0.2353 1.0000 0.5765 0.7569 0.0000 0.9804 0.8000 1.0000 1.0000 0.7608 0.8000 0.0667 0.6314 0.7569 1.0000 0.8275 0.1804 0.2000 0.8510 0.4706 0.3608 0.8353 0.1059 1.0000 1.0000 0.1608 0.3804 1.0000 0.5059 0.7412 0.1922 0.7451 0.7373 1.0000 0.8471 0.1608 0.5451 0.7412 0.8157 0.7059 0.6314 0.0863 1.0000 0.7765 0.1804 0.6902 0.2118 0.6510 0.6275 0.2824 0.8902 0.3922 0.2902 0.4902 0.3255 0.6863 0.0510 0.4510 0.4157 0.5608 0.6706 0.0510 0.4471 0.3804 0.4706 0.3804 0.3059 0.4314 0.4275 0.5804 0.6157 0.2902 0.2667 0.3569 0.2510 0.1451 0.2157 0.1255 0.2157 0.3608 0.2314 0.1608 0.0510 0.2157 0.0157 0.1412 0.2314 0.2157 0.1608 0.2706 0.0706 0.0157 0.2353 0.1451 0.1608 0.3059 0.0000 0.0000 0.0353 0.0471 0.1255 0.0902 0.0471 0.0157 0.0353 0.0902 0.0314 0.0118 0.0000 0.0157 0.0000 0.0000 0.0118 0.0157 0.0157 0.0314 0.0000 0.0706 0.0314 0.0000 0.0510 0.0353 0.0667 0.1608 0.2157 0.1020 0.1608 0.2510 0.2314 0.1373 0.1059 0.2510 0.3059 0.2314 0.2902 0.0902 0.2000 0.1255 0.1608 0.3804 0.2118 0.3020 0.2510 0.1255 0.2510 0.0510 0.3216 0.4157 0.3255 0.4471 0.3255 0.1255 0.0706 0.2353 0.0706 0.0902 0.3216 0.1451 0.1412 0.2510 0.5451 0.6706 0.0706 0.0353 0.7059 0.1059 0.6000 0.0706 0.0275 0.8000 0.6902 0.5804 0.7216 0.4902 0.4314 0.0157 0.1373 0.2157 0.1059 0.7922 0.0510 0.6314 0.7059 0.5059 0.2000 0.8000 0.4667 0.0667 0.3922 0.6902 0.2510 0.5059 0.8510 0.6353 0.7569 0.8824 1.0000 0.4510 0.2510 0.1255 0.0118 0.4000 0.1020 1.0000 0.3804 1.0000 1.0000 0.0157 1.0000 0.2706 1.0000 1.0000 1.0000 0.5216 0.9059 0.7608 1.0000 1.0000 0.3765 1.0000 1.0000 1.0000 0.7804 0.2706 0.4275 1.0000 0.3255 0.8471 0.6157 1.0000 1.0000 1.0000 0.0353 1.0000 0.7059 1.0000 1.0000 1.0000 0.4510];
    case 'random - constant'
        red=[ 0.2353 1.0000 1.0000 1.0000 0.9412 0.9608 1.0000 0.5608 0.9059 0.1451 0.9922 1.0000 1.0000 0.8157 0.8824 1.0000 1.0000 0.9608 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.3451 1.0000 1.0000 0.6275 1.0000 0.7059 0.6353 0.7255 1.0000 1.0000 0.9569 1.0000 1.0000 0.2510 1.0000 0.7608 1.0000 1.0000 1.0000 1.0000 0.1451 0.8000 1.0000 1.0000 0.3608 0.4157 1.0000 1.0000 1.0000 0.7608 0.6314 1.0000 1.0000 1.0000 1.0000 0.3255 1.0000 0.5922 0.7765 0.1608 0.2902 0.4471 1.0000 0.1059 1.0000 0.7608 0.8353 0.8902 0.2471 0.2706 1.0000 0.9608 0.4353 0.5922 1.0000 1.0000 0.1804 1.0000 1.0000 0.9804 1.0000 0.7608 1.0000 0.4863 1.0000 1.0000 1.0000 0.9608 0.2902 0.8314 0.8706 0.3020 0.9608 0.8706 0.7059 0.3804 0.2314 1.0000 0.7412 0.0902 1.0000 1.0000 1.0000 1.0000 0.0157 0.3608 0.5804 1.0000 0.3059 0.9059 0.5804 1.0000 0.3608 1.0000 0.5255 0.8706 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.5059 1.0000 0.7765 1.0000 0.2706 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.0157 0.7255 1.0000 1.0000 0.7373 0.8706 1.0000 1.0000 0.8353 1.0000 0.6510 1.0000 1.0000 0.1804 1.0000 0.4157 1.0000 0.8000 0.7412 1.0000 0.1020 0.6118 1.0000 1.0000 1.0000 1.0000 0.1255 1.0000 0.5608 0.8902 0.3961 1.0000 1.0000 0.2902 1.0000 1.0000 0.6902 0.5804 0.9961 1.0000 1.0000 1.0000 1.0000 0.9412 1.0000 0.2275 0.1412 1.0000 1.0000 1.0000 0.2000 0.2902 1.0000 0.2118 1.0000 0.9059 0.4863 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.3059 1.0000 1.0000 0.0706 1.0000 1.0000 1.0000 1.0000 1.0000 0.5804 0.7961 1.0000 0.5059 1.0000 1.0000 0.6000 0.5804 0.7059 1.0000 0.4000 1.0000 0.9608 1.0000 1.0000 1.0000 0.3059 0.7216 0.3255 1.0000 1.0000 1.0000 1.0000 0.6000 0.2902 0.3608 0.6157 1.0000 0.4157 0.4863 1.0000 0.9608 1.0000 1.0000 1.0000 0.6353 0.1451 0.2000 0.5373 0.9451 0.0157 0.9451 0.6157 0.7804];
        green=[ 1.0000 0.0353 0.3608 0.4510 1.0000 1.0000 1.0000 1.0000 0.8510 1.0000 0.3608 0.4314 0.2902 1.0000 0.6706 0.9804 0.3216 1.0000 0.9059 0.7804 1.0000 1.0000 0.7608 0.7804 0.0000 1.0000 1.0000 1.0000 0.5059 0.7059 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.3412 0.9059 0.7608 0.5412 1.0000 1.0000 0.3059 0.6157 1.0000 1.0000 1.0000 0.7922 0.9216 1.0000 0.0510 0.1804 0.6000 1.0000 1.0000 0.0510 0.9412 0.9451 0.6353 1.0000 1.0000 0.6706 0.9922 1.0000 0.7451 1.0000 0.7451 1.0000 1.0000 1.0000 1.0000 0.9608 1.0000 1.0000 1.0000 0.3608 1.0000 0.8471 0.1059 0.3961 1.0000 0.3451 1.0000 0.8471 0.3569 0.9804 1.0000 1.0000 0.6706 0.2353 0.4510 1.0000 1.0000 1.0000 1.0000 0.9608 0.8824 0.8863 0.5255 1.0000 1.0000 1.0000 1.0000 1.0000 0.8824 0.3451 1.0000 1.0000 0.3255 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.6353 1.0000 0.0157 0.3608 1.0000 0.4157 1.0000 1.0000 0.8157 0.8314 0.9216 0.7608 1.0000 0.4706 0.8510 1.0000 1.0000 0.4000 0.9255 0.0706 0.9059 0.7059 0.8510 0.2706 0.4471 0.7059 0.4118 0.6118 1.0000 1.0000 0.3059 0.8706 0.7608 0.6353 0.6902 0.2510 0.7804 1.0000 0.1961 0.8000 0.2157 1.0000 0.0706 0.8902 1.0000 1.0000 1.0000 1.0000 0.0000 0.9804 0.9922 1.0000 0.9608 0.6667 0.7608 0.1216 0.8902 1.0000 1.0000 0.4118 1.0000 1.0000 1.0000 0.7412 0.1765 0.7608 1.0000 0.3804 1.0000 1.0000 0.4902 0.1059 0.6902 1.0000 1.0000 1.0000 0.7608 1.0000 0.8157 0.3255 1.0000 0.4902 0.1255 1.0000 0.2000 0.7059 0.9059 1.0000 0.8510 1.0000 1.0000 0.9922 0.6902 0.3412 0.8902 0.3608 0.2353 0.4510 0.5373 0.7255 0.6157 0.2706 1.0000 1.0000 1.0000 0.6863 1.0000 0.9569 0.4157 1.0000 1.0000 0.8314 1.0000 1.0000 1.0000 1.0000 1.0000 0.5804 0.9804 0.6510 1.0000 1.0000 1.0000 0.0824 1.0000 1.0000 0.9059 1.0000 0.6353 1.0000 0.9961 1.0000 1.0000 1.0000 0.6000 1.0000 0.0471 1.0000 1.0000 1.0000];
        blue=[ 1.0000 1.0000 0.2275 1.0000 0.8157 0.8353 0.7451 0.5804 1.0000 1.0000 1.0000 1.0000 1.0000 0.8510 1.0000 0.8157 0.4157 0.9608 0.8353 1.0000 0.5922 0.8314 0.8314 0.9961 1.0000 0.5765 0.0863 0.5451 1.0000 0.6314 1.0000 0.8353 0.9451 0.0706 0.0157 0.7608 0.0902 0.9216 1.0000 0.9255 1.0000 0.5373 0.5020 1.0000 1.0000 1.0000 0.4510 0.1020 1.0000 1.0000 0.8510 1.0000 0.5608 1.0000 0.5804 0.9608 1.0000 0.0706 0.2706 1.0000 0.7608 0.3255 1.0000 1.0000 1.0000 1.0000 0.8157 0.2863 0.5451 0.7255 0.7608 0.8510 1.0000 0.7608 1.0000 0.1961 1.0000 0.8510 1.0000 0.5804 0.7373 1.0000 1.0000 0.7020 1.0000 1.0000 1.0000 0.5804 0.9922 1.0000 1.0000 0.4275 0.6863 1.0000 1.0000 0.7608 1.0000 1.0000 1.0000 1.0000 0.5059 1.0000 0.6471 0.9608 1.0000 0.8706 1.0000 0.0824 0.6157 1.0000 1.0000 0.8000 0.2824 0.8000 1.0000 1.0000 1.0000 1.0000 0.7373 1.0000 0.8000 1.0000 0.1922 0.6863 1.0000 0.0824 0.9608 0.8157 1.0000 1.0000 1.0000 0.4706 0.9059 1.0000 0.8510 0.9451 0.2863 0.6706 1.0000 0.9608 1.0000 1.0000 1.0000 1.0000 0.0902 0.7961 0.8314 1.0000 1.0000 0.9922 1.0000 1.0000 0.7059 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.7059 1.0000 0.1255 0.4353 0.1804 0.8824 0.1451 1.0000 0.8314 1.0000 1.0000 1.0000 0.0314 0.7608 1.0000 0.3373 1.0000 0.9608 1.0000 0.8706 1.0000 1.0000 1.0000 0.7059 1.0000 0.8510 0.9804 1.0000 1.0000 0.0902 0.2000 0.6824 0.7608 0.9216 1.0000 0.3059 1.0000 0.8706 0.8667 1.0000 0.4000 1.0000 0.5569 0.9451 1.0000 0.9412 0.1059 1.0000 0.1216 0.2275 0.3255 0.8667 1.0000 1.0000 1.0000 0.3059 1.0000 0.7569 1.0000 0.9451 0.6824 0.8275 0.2510 0.5216 0.8902 1.0000 0.4353 0.5373 0.9922 1.0000 1.0000 1.0000 0.3608 0.8353 0.2824 0.9412 1.0000 0.9451 1.0000 1.0000 0.9804 1.0000 0.8902 1.0000 0.4863 1.0000 0.9804 0.2353 1.0000 1.0000 1.0000 1.0000 0.8353 1.0000 0.8353 1.0000 0.9451];
    case 'random - linear'
        red=[ 0.0000 0.0000 0.0000 0.0000 0.0157 0.0157 0.0157 0.0353 0.0353 0.0706 0.0353 0.0275 0.0157 0.0706 0.0510 0.0275 0.0000 0.0706 0.0000 0.0706 0.0706 0.0510 0.0510 0.0314 0.0000 0.1451 0.0353 0.0275 0.1412 0.0157 0.1608 0.1608 0.2000 0.1059 0.0510 0.1451 0.1255 0.0863 0.2118 0.1059 0.2000 0.1059 0.0510 0.1451 0.1412 0.3059 0.1804 0.1804 0.1804 0.2157 0.2824 0.1255 0.1059 0.1451 0.2353 0.3804 0.1765 0.0353 0.0000 0.1765 0.3804 0.1373 0.3608 0.3255 0.4314 0.3804 0.3451 0.1804 0.3804 0.2471 0.3608 0.3804 0.3608 0.4157 0.3569 0.1059 0.3020 0.4314 0.3451 0.0157 0.0000 0.7059 0.1059 0.2510 0.3373 0.2353 0.5804 0.3255 0.5255 0.1804 0.3255 0.0510 0.3608 0.9059 0.3961 0.4000 0.4902 0.4000 0.5059 0.4000 0.5451 0.5804 0.0706 0.6824 0.7059 0.3255 0.2667 0.3020 0.4157 0.6157 0.8353 0.5255 0.3608 0.6000 0.5255 0.5216 0.3608 0.9804 0.2157 0.4824 0.4824 0.4157 0.2706 0.4314 0.3961 0.2902 0.3765 0.2902 0.9020 0.2353 0.7451 0.3608 0.9451 0.3804 0.3059 0.0471 0.0902 0.2157 0.3922 0.3451 0.8353 0.5804 0.4157 0.4118 0.6510 0.7255 0.3608 0.5020 0.7059 0.1451 0.7608 0.3451 0.0275 0.8510 0.0510 1.0000 0.4863 1.0000 0.7255 0.1451 1.0000 0.7804 0.4000 0.5020 0.0157 0.3608 0.9961 0.6000 0.8510 0.7608 0.8000 0.0510 0.5059 1.0000 0.6157 0.3412 0.8706 1.0000 0.6902 0.6824 0.0353 0.6902 0.4863 0.7412 0.2118 0.8157 1.0000 0.1608 0.4157 0.0353 0.9961 0.9059 0.3608 1.0000 0.0314 0.7922 1.0000 0.2824 0.5922 0.2706 0.4902 0.0275 0.5804 1.0000 0.2824 0.6275 1.0000 0.5373 0.2275 0.3059 0.2902 0.4353 0.9765 0.9451 0.4902 1.0000 0.7059 0.5216 1.0000 1.0000 1.0000 0.6510 1.0000 0.7569 0.8902 0.5255 0.7569 0.5255 1.0000 1.0000 1.0000 0.7608 0.8000 0.1451 0.7059 1.0000 1.0000 1.0000 1.0000 0.7765 1.0000 1.0000 0.0902 0.9608 0.8157 0.8902 0.2824 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
        green=[ 0.0000 0.0000 0.0000 0.0157 0.0118 0.0118 0.0157 0.0000 0.0471 0.0000 0.0706 0.0510 0.0902 0.0000 0.0706 0.0510 0.0902 0.0510 0.0902 0.0706 0.0510 0.0157 0.1020 0.1451 0.2000 0.0000 0.0000 0.0902 0.1451 0.1451 0.0902 0.0510 0.0157 0.1020 0.0902 0.0353 0.0157 0.1922 0.1451 0.1804 0.2471 0.1451 0.1255 0.2471 0.3412 0.1451 0.1373 0.1373 0.3216 0.1961 0.0863 0.3373 0.2510 0.3451 0.1059 0.0275 0.2902 0.2314 0.2353 0.2902 0.0275 0.1804 0.3451 0.2471 0.2353 0.3020 0.1373 0.2706 0.0863 0.2275 0.0706 0.0824 0.3059 0.1059 0.2510 0.1255 0.4157 0.1373 0.3216 0.5020 0.5255 0.2000 0.6510 0.2902 0.4353 0.5216 0.3451 0.3216 0.1451 0.5412 0.6902 0.4902 0.3059 0.0667 0.3412 0.2510 0.3804 0.4667 0.4902 0.4157 0.1059 0.3255 0.3255 0.0863 0.3569 0.4471 0.8000 0.3608 0.2863 0.5608 0.3804 0.2902 0.2157 0.2157 0.4275 0.3804 0.7216 0.3451 0.6510 0.4902 0.4314 0.5804 0.1451 0.4510 0.8000 0.5059 0.5608 0.6000 0.4902 0.9765 0.6510 0.0510 0.0353 0.7608 0.5922 0.9804 0.5765 0.6824 0.7255 0.7255 0.7059 0.5922 0.7255 1.0000 0.1255 0.1451 0.7412 0.9059 0.7804 1.0000 0.7412 1.0000 0.8353 0.3255 1.0000 0.7373 0.8157 0.5765 1.0000 0.7451 0.0157 0.1255 0.0824 0.5608 1.0000 0.6510 0.6510 0.1804 0.6706 1.0000 0.7216 0.9608 0.7059 0.2863 0.3412 1.0000 0.4863 0.4353 0.6510 1.0000 1.0000 0.7451 0.6157 0.9922 0.6824 0.6157 0.9255 1.0000 0.8157 0.2824 0.3804 0.5255 1.0000 0.5412 0.9059 1.0000 0.3451 1.0000 0.9608 0.5804 1.0000 1.0000 0.9255 0.5804 1.0000 0.0000 0.2157 0.8118 0.9804 1.0000 1.0000 1.0000 1.0000 1.0000 0.9765 0.9569 0.9255 1.0000 0.1255 0.3922 0.1451 0.9255 0.1255 0.9059 1.0000 0.4353 0.5216 1.0000 0.6157 0.5451 0.7804 0.8902 0.7608 1.0000 0.9608 1.0000 0.1765 0.2000 0.4471 1.0000 0.4510 0.2863 1.0000 0.7216 1.0000 0.8667 0.9608 0.2314 0.4118 0.2902 1.0000 0.3804 1.0000 0.9059 0.3373 0.5804];
        blue=[ 0.0000 0.0000 0.0000 0.0000 0.0157 0.0157 0.0157 0.0353 0.0000 0.0157 0.0000 0.0353 0.0157 0.0667 0.0157 0.0863 0.0902 0.0706 0.1059 0.0706 0.0902 0.1608 0.0902 0.0902 0.0510 0.1255 0.2510 0.1804 0.0275 0.1608 0.0902 0.1373 0.1412 0.1608 0.2353 0.2157 0.2667 0.1451 0.0706 0.1608 0.0157 0.2157 0.2902 0.1059 0.0275 0.0706 0.2000 0.2275 0.0510 0.1451 0.2157 0.1255 0.2314 0.1255 0.2706 0.2157 0.1804 0.3961 0.4353 0.2157 0.2863 0.3804 0.0118 0.1608 0.0706 0.0706 0.2863 0.3059 0.3255 0.3255 0.3608 0.3608 0.1608 0.3255 0.2510 0.6314 0.1608 0.3255 0.2353 0.4000 0.4000 0.0275 0.1961 0.4157 0.2118 0.2353 0.0824 0.3608 0.3451 0.3059 0.0275 0.5059 0.4000 0.1020 0.3608 0.4353 0.2510 0.2706 0.1373 0.3255 0.5059 0.2706 0.7804 0.4353 0.1451 0.4510 0.1765 0.5765 0.5569 0.0902 0.0510 0.4706 0.7255 0.4902 0.3765 0.4353 0.2706 0.0275 0.5059 0.4157 0.4902 0.4157 1.0000 0.5451 0.2510 0.6510 0.5255 0.5804 0.1020 0.2902 0.1059 1.0000 0.5608 0.4000 0.6510 0.5451 0.9255 0.6902 0.4902 0.5451 0.0902 0.4706 0.5216 0.2471 0.9020 0.8157 0.6118 0.3059 0.2353 0.5922 0.2510 0.0314 0.9059 0.5961 0.1412 0.0824 0.5255 0.2471 0.1059 0.9608 0.5451 0.9608 1.0000 0.8353 0.7059 0.9059 0.3020 1.0000 0.4353 0.2000 0.4667 0.9922 0.7804 0.3608 1.0000 0.2157 0.7059 0.0471 0.7373 0.2902 0.2706 0.6706 1.0000 0.4157 1.0000 0.7255 0.1922 0.6902 0.9608 1.0000 0.8353 0.8000 0.8471 0.2000 1.0000 0.4706 0.8471 0.8706 0.7608 1.0000 0.6667 1.0000 0.8510 0.1569 0.9412 1.0000 0.7059 1.0000 1.0000 1.0000 1.0000 0.6510 0.4000 0.4157 1.0000 0.5059 0.8902 0.6667 0.9412 1.0000 1.0000 1.0000 1.0000 0.9569 0.6824 1.0000 1.0000 0.9020 0.1569 0.8706 0.5765 1.0000 1.0000 1.0000 1.0000 0.2353 0.9608 0.4706 0.3059 0.9373 0.6157 1.0000 0.5922 1.0000 0.4510 1.0000 1.0000 0.7059 0.8118 0.8471 0.5059 1.0000 0.1373 1.0000 0.7804 1.0000];
    case 'bipolar 5'
        red=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0353 0.0471 0.0510 0.0706 0.0902 0.1059 0.1216 0.1255 0.1451 0.1608 0.1804 0.1922 0.2000 0.2157 0.2353 0.2510 0.2706 0.2863 0.3020 0.3059 0.3255 0.3451 0.3608 0.3765 0.3922 0.4000 0.4157 0.4353 0.4510 0.4706 0.4863 0.5020 0.5059 0.5373 0.5451 0.5608 0.5804 0.6000 0.6157 0.6314 0.6471 0.6667 0.6824 0.6902 0.7059 0.7255 0.7451 0.7608 0.7804 0.7961 0.8118 0.8314 0.8471 0.8510 0.8706 0.8902 0.9059 0.9255 0.9451 0.9608 0.9804 1.0000 0.9922 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9804 0.9765 0.9608 0.9451 0.9255 0.9059 0.9020 0.8902 0.8706 0.8510 0.8471 0.8314 0.8157 0.8000 0.7804 0.7765 0.7608 0.7451 0.7255 0.7059 0.7020 0.6902 0.6706 0.6510 0.6471 0.6314 0.6157 0.6000 0.5804];
        green=[ 0.3451 0.3608 0.3608 0.3804 0.3804 0.3961 0.4000 0.4118 0.4157 0.4314 0.4353 0.4471 0.4510 0.4510 0.4706 0.4824 0.4902 0.4902 0.5059 0.5059 0.5216 0.5255 0.5412 0.5451 0.5569 0.5608 0.5608 0.5804 0.5922 0.6000 0.6000 0.6157 0.6157 0.6314 0.6353 0.6510 0.6510 0.6667 0.6706 0.6824 0.6902 0.7020 0.7059 0.7059 0.7255 0.7255 0.7412 0.7451 0.7608 0.7608 0.7765 0.7804 0.7922 0.8000 0.8000 0.8157 0.8275 0.8353 0.8353 0.8510 0.8510 0.8667 0.8706 0.8863 0.8902 0.8902 0.8902 0.8902 0.9020 0.9020 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9216 0.9216 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9373 0.9373 0.9373 0.9412 0.9412 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9569 0.9569 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9765 0.9765 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9922 0.9922 0.9922 0.9961 0.9961 0.9961 1.0000 0.9922 0.9608 0.9608 0.9451 0.9373 0.9216 0.9059 0.8902 0.8863 0.8706 0.8510 0.8471 0.8314 0.8157 0.8000 0.7922 0.7804 0.7608 0.7451 0.7373 0.7255 0.7059 0.6902 0.6824 0.6706 0.6510 0.6353 0.6275 0.6118 0.6000 0.5804 0.5608 0.5451 0.5412 0.5255 0.5059 0.4902 0.4706 0.4510 0.4471 0.4314 0.4157 0.3961 0.3804 0.3608 0.3451 0.3255 0.3059 0.2902 0.2706 0.2510 0.2353 0.2157 0.2000 0.1804 0.1608 0.1412 0.1216 0.0902 0.0706 0.0510 0.0314 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000];
        blue=[ 0.3451 0.3608 0.3608 0.3804 0.3804 0.3961 0.4000 0.4118 0.4157 0.4314 0.4353 0.4471 0.4510 0.4510 0.4706 0.4824 0.4902 0.4902 0.5059 0.5059 0.5216 0.5255 0.5412 0.5451 0.5569 0.5608 0.5608 0.5804 0.5922 0.6000 0.6000 0.6157 0.6157 0.6314 0.6353 0.6510 0.6510 0.6667 0.6706 0.6824 0.6902 0.7020 0.7059 0.7059 0.7255 0.7255 0.7412 0.7451 0.7608 0.7608 0.7765 0.7804 0.7922 0.8000 0.8000 0.8157 0.8275 0.8353 0.8353 0.8510 0.8510 0.8667 0.8706 0.8863 0.8902 0.8902 0.8902 0.8902 0.9020 0.9020 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9216 0.9216 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9373 0.9373 0.9373 0.9412 0.9412 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9569 0.9569 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9765 0.9765 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9922 0.9922 0.9922 0.9961 0.9961 0.9961 1.0000 0.9922 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9804 0.9765 0.9608 0.9451 0.9255 0.9059 0.9020 0.8902 0.8706 0.8510 0.8471 0.8314 0.8157 0.8000 0.7804 0.7765 0.7608 0.7451 0.7255 0.7059 0.7020 0.6902 0.6706 0.6510 0.6471 0.6314 0.6157 0.6000 0.5804];
    case 'bipolar 4'
        red=[ 0.2471 0.2510 0.2706 0.2902 0.3059 0.3255 0.3412 0.3451 0.3608 0.3804 0.4000 0.4118 0.4275 0.4353 0.4510 0.4706 0.4824 0.4902 0.5059 0.5255 0.5373 0.5451 0.5608 0.5765 0.5804 0.6000 0.6157 0.6275 0.6353 0.6510 0.6510 0.6706 0.6863 0.6902 0.7059 0.7216 0.7255 0.7451 0.7569 0.7608 0.7765 0.7804 0.8000 0.8000 0.8157 0.8275 0.8353 0.8510 0.8510 0.8706 0.8706 0.8863 0.8902 0.9059 0.9059 0.9255 0.9255 0.9412 0.9451 0.9608 0.9608 0.9765 0.9804 0.9922 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9922 0.9765 0.9608 0.9451 0.9255 0.9059 0.8902 0.8824 0.8667 0.8510 0.8353 0.8157 0.8000 0.7804 0.7608 0.7569 0.7412 0.7255 0.7059 0.6902 0.6706 0.6510 0.6353 0.6314 0.6157 0.6000 0.5804 0.5608 0.5451 0.5255 0.5059 0.5059 0.4902 0.4706 0.4510 0.4353 0.4157 0.4000 0.3922 0.3804 0.3608 0.3451 0.3255 0.3059 0.2902 0.2824 0.2667 0.2510 0.2353 0.2157 0.2000 0.1922 0.1608 0.1569 0.1412 0.1255 0.1059 0.0902 0.0824 0.0667 0.0471 0.0314 0.0157 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000];
        green=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0353 0.0510 0.0510 0.0706 0.0902 0.1059 0.1255 0.1412 0.1608 0.1765 0.1804 0.2000 0.2157 0.2353 0.2510 0.2667 0.2863 0.3020 0.3059 0.3255 0.3451 0.3608 0.3804 0.3922 0.4118 0.4275 0.4353 0.4510 0.4706 0.4902 0.5059 0.5059 0.5373 0.5451 0.5608 0.5804 0.6000 0.6157 0.6314 0.6353 0.6510 0.6706 0.6902 0.7059 0.7255 0.7412 0.7569 0.7608 0.7804 0.8000 0.8157 0.8314 0.8510 0.8667 0.8824 0.8902 0.9059 0.9255 0.9451 0.9569 0.9765 0.9922 0.9922 0.9765 0.9608 0.9451 0.9255 0.9059 0.8902 0.8824 0.8667 0.8510 0.8353 0.8157 0.8000 0.7804 0.7608 0.7569 0.7412 0.7255 0.7059 0.6902 0.6706 0.6510 0.6353 0.6314 0.6157 0.6000 0.5804 0.5608 0.5451 0.5255 0.5059 0.5059 0.4902 0.4706 0.4510 0.4353 0.4157 0.4000 0.3922 0.3804 0.3608 0.3451 0.3255 0.3059 0.2902 0.2824 0.2667 0.2510 0.2353 0.2157 0.2000 0.1922 0.1608 0.1569 0.1412 0.1255 0.1059 0.0902 0.0824 0.0667 0.0471 0.0314 0.0157 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000];
        blue=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0353 0.0510 0.0510 0.0706 0.0902 0.1059 0.1255 0.1412 0.1608 0.1765 0.1804 0.2000 0.2157 0.2353 0.2510 0.2667 0.2863 0.3020 0.3059 0.3255 0.3451 0.3608 0.3804 0.3922 0.4118 0.4275 0.4353 0.4510 0.4706 0.4902 0.5059 0.5059 0.5373 0.5451 0.5608 0.5804 0.6000 0.6157 0.6314 0.6353 0.6510 0.6706 0.6902 0.7059 0.7255 0.7412 0.7569 0.7608 0.7804 0.8000 0.8157 0.8314 0.8510 0.8667 0.8824 0.8902 0.9059 0.9255 0.9451 0.9569 0.9765 0.9922 0.9922 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9922 0.9804 0.9765 0.9608 0.9608 0.9451 0.9412 0.9255 0.9255 0.9059 0.9059 0.8902 0.8863 0.8706 0.8706 0.8510 0.8510 0.8353 0.8275 0.8157 0.8000 0.8000 0.7804 0.7765 0.7608 0.7569 0.7451 0.7255 0.7216 0.7059 0.6902 0.6863 0.6706 0.6510 0.6510 0.6353 0.6275 0.6157 0.6000 0.5804 0.5765 0.5608 0.5451 0.5373 0.5255 0.5059 0.4902 0.4824 0.4706 0.4510 0.4353 0.4275 0.4118 0.4000 0.3804 0.3608 0.3451 0.3412 0.3255 0.3059 0.2902 0.2706 0.2510 0.2471];
    case 'bipolar 3'
        red=[ 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9804 0.9608 0.9608 0.9569 0.9451 0.9412 0.9255 0.9255 0.9059 0.9059 0.9020 0.8902 0.8863 0.8706 0.8706 0.8510 0.8510 0.8471 0.8353 0.8314 0.8157 0.8157 0.8000 0.8000 0.7922 0.7804 0.7765 0.7608 0.7608 0.7451 0.7451 0.7373 0.7255 0.7216 0.7059 0.7059 0.6902 0.6902 0.6824 0.6706 0.6667 0.6510 0.6510 0.6353 0.6353 0.6275 0.6157 0.6118 0.6000 0.5961 0.5804 0.5804 0.5608 0.5608 0.5569 0.5451 0.5412 0.5255 0.5255 0.5059 0.5059 0.5020 0.4902 0.4863 0.4706 0.4706 0.4510 0.4510 0.4471 0.4353 0.4314 0.4157 0.4157 0.4000 0.4000 0.3922 0.3804 0.3765 0.3608 0.3608 0.3451 0.3451 0.3373 0.3255 0.3216 0.3059 0.3059 0.2902 0.2902 0.2824 0.2706 0.2667 0.2510 0.2510 0.2353 0.2353 0.2275 0.2157 0.2118 0.2000 0.1961 0.1804 0.1804 0.1608 0.1608 0.1569 0.1451 0.1412 0.1255 0.1255 0.1059 0.1059 0.1020 0.0902 0.0863 0.0706 0.0706 0.0510 0.0510 0.0471 0.0353 0.0314 0.0157 0.0157 0.0000 0.0000 0.0000];
        green=[ 0.0000 0.0000 0.0157 0.0157 0.0314 0.0353 0.0471 0.0510 0.0510 0.0706 0.0706 0.0863 0.0902 0.1020 0.1059 0.1059 0.1255 0.1255 0.1412 0.1451 0.1569 0.1608 0.1608 0.1804 0.1804 0.1961 0.2000 0.2118 0.2157 0.2275 0.2353 0.2353 0.2510 0.2510 0.2667 0.2706 0.2824 0.2902 0.2902 0.3059 0.3059 0.3216 0.3255 0.3373 0.3451 0.3451 0.3608 0.3608 0.3765 0.3804 0.3922 0.4000 0.4000 0.4157 0.4157 0.4314 0.4353 0.4471 0.4510 0.4510 0.4706 0.4706 0.4863 0.4902 0.5020 0.5059 0.5059 0.5255 0.5255 0.5412 0.5451 0.5569 0.5608 0.5608 0.5804 0.5804 0.5961 0.6000 0.6118 0.6157 0.6275 0.6353 0.6353 0.6510 0.6510 0.6667 0.6706 0.6824 0.6902 0.6902 0.7059 0.7059 0.7216 0.7255 0.7373 0.7451 0.7451 0.7608 0.7608 0.7765 0.7804 0.7922 0.8000 0.8000 0.8157 0.8157 0.8314 0.8353 0.8471 0.8510 0.8510 0.8706 0.8706 0.8863 0.8902 0.9020 0.9059 0.9059 0.9255 0.9255 0.9412 0.9451 0.9569 0.9608 0.9608 0.9804 0.9804 1.0000 0.9804 0.9804 0.9608 0.9608 0.9569 0.9451 0.9412 0.9255 0.9255 0.9059 0.9059 0.9020 0.8902 0.8863 0.8706 0.8706 0.8510 0.8510 0.8471 0.8353 0.8314 0.8157 0.8157 0.8000 0.8000 0.7922 0.7804 0.7765 0.7608 0.7608 0.7451 0.7451 0.7373 0.7255 0.7216 0.7059 0.7059 0.6902 0.6902 0.6824 0.6706 0.6667 0.6510 0.6510 0.6353 0.6353 0.6275 0.6157 0.6118 0.6000 0.5961 0.5804 0.5804 0.5608 0.5608 0.5569 0.5451 0.5412 0.5255 0.5255 0.5059 0.5059 0.5020 0.4902 0.4863 0.4706 0.4706 0.4510 0.4510 0.4471 0.4353 0.4314 0.4157 0.4157 0.4000 0.4000 0.3922 0.3804 0.3765 0.3608 0.3608 0.3451 0.3451 0.3373 0.3255 0.3216 0.3059 0.3059 0.2902 0.2902 0.2824 0.2706 0.2667 0.2510 0.2510 0.2353 0.2353 0.2275 0.2157 0.2118 0.2000 0.1961 0.1804 0.1804 0.1608 0.1608 0.1569 0.1451 0.1412 0.1255 0.1255 0.1059 0.1059 0.1020 0.0902 0.0863 0.0706 0.0706 0.0510 0.0510 0.0471 0.0353 0.0314 0.0157 0.0157 0.0000 0.0000 0.0000];
        blue=[ 0.0000 0.0000 0.0157 0.0157 0.0314 0.0353 0.0471 0.0510 0.0510 0.0706 0.0706 0.0863 0.0902 0.1020 0.1059 0.1059 0.1255 0.1255 0.1412 0.1451 0.1569 0.1608 0.1608 0.1804 0.1804 0.1961 0.2000 0.2118 0.2157 0.2275 0.2353 0.2353 0.2510 0.2510 0.2667 0.2706 0.2824 0.2902 0.2902 0.3059 0.3059 0.3216 0.3255 0.3373 0.3451 0.3451 0.3608 0.3608 0.3765 0.3804 0.3922 0.4000 0.4000 0.4157 0.4157 0.4314 0.4353 0.4471 0.4510 0.4510 0.4706 0.4706 0.4863 0.4902 0.5020 0.5059 0.5059 0.5255 0.5255 0.5412 0.5451 0.5569 0.5608 0.5608 0.5804 0.5804 0.5961 0.6000 0.6118 0.6157 0.6275 0.6353 0.6353 0.6510 0.6510 0.6667 0.6706 0.6824 0.6902 0.6902 0.7059 0.7059 0.7216 0.7255 0.7373 0.7451 0.7451 0.7608 0.7608 0.7765 0.7804 0.7922 0.8000 0.8000 0.8157 0.8157 0.8314 0.8353 0.8471 0.8510 0.8510 0.8706 0.8706 0.8863 0.8902 0.9020 0.9059 0.9059 0.9255 0.9255 0.9412 0.9451 0.9569 0.9608 0.9608 0.9804 0.9804 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
    case 'bipolar 2'
        red=[ 0.5020 0.5059 0.5059 0.5255 0.5255 0.5412 0.5451 0.5569 0.5608 0.5608 0.5765 0.5804 0.5961 0.6000 0.6118 0.6157 0.6275 0.6353 0.6353 0.6510 0.6510 0.6667 0.6706 0.6824 0.6902 0.6902 0.7020 0.7059 0.7216 0.7255 0.7373 0.7451 0.7451 0.7608 0.7608 0.7765 0.7804 0.7922 0.8000 0.8000 0.8157 0.8157 0.8275 0.8353 0.8471 0.8510 0.8510 0.8706 0.8706 0.8863 0.8902 0.9020 0.9059 0.9059 0.9255 0.9255 0.9412 0.9451 0.9569 0.9608 0.9608 0.9804 0.9804 0.9961 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9804 0.9608 0.9451 0.9373 0.9216 0.9059 0.8902 0.8706 0.8510 0.8353 0.8275 0.8118 0.7961 0.7804 0.7608 0.7451 0.7255 0.7059 0.7020 0.6863 0.6706 0.6510 0.6353 0.6157 0.6000 0.5922 0.5765 0.5608 0.5451 0.5255 0.5059 0.4902 0.4706 0.4667 0.4471 0.4353 0.4157 0.4000 0.3804 0.3608 0.3451 0.3412 0.3216 0.3059 0.2902 0.2706 0.2510 0.2471 0.2275 0.2157 0.1961 0.1804 0.1608 0.1451 0.1255 0.1216 0.1020 0.0902 0.0706 0.0510 0.0353 0.0275 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000];
        green=[ 0.0000 0.0157 0.0314 0.0471 0.0510 0.0706 0.0902 0.1059 0.1255 0.1412 0.1569 0.1608 0.1804 0.2000 0.2157 0.2353 0.2510 0.2667 0.2824 0.2902 0.3059 0.3255 0.3451 0.3608 0.3765 0.3922 0.4000 0.4157 0.4353 0.4510 0.4706 0.4863 0.5020 0.5059 0.5255 0.5451 0.5608 0.5765 0.5961 0.6118 0.6275 0.6353 0.6510 0.6706 0.6902 0.7020 0.7216 0.7373 0.7451 0.7608 0.7804 0.8000 0.8157 0.8275 0.8471 0.8510 0.8706 0.8902 0.9059 0.9255 0.9412 0.9569 0.9608 0.9804 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9804 0.9608 0.9451 0.9373 0.9216 0.9059 0.8902 0.8706 0.8510 0.8353 0.8275 0.8118 0.7961 0.7804 0.7608 0.7451 0.7255 0.7059 0.7020 0.6863 0.6706 0.6510 0.6353 0.6157 0.6000 0.5922 0.5765 0.5608 0.5451 0.5255 0.5059 0.4902 0.4706 0.4667 0.4471 0.4353 0.4157 0.4000 0.3804 0.3608 0.3451 0.3412 0.3216 0.3059 0.2902 0.2706 0.2510 0.2471 0.2275 0.2157 0.1961 0.1804 0.1608 0.1451 0.1255 0.1216 0.1020 0.0902 0.0706 0.0510 0.0353 0.0275 0.0000];
        blue=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0510 0.0706 0.0902 0.1059 0.1255 0.1412 0.1569 0.1608 0.1804 0.2000 0.2157 0.2353 0.2510 0.2667 0.2824 0.2902 0.3059 0.3255 0.3451 0.3608 0.3765 0.3922 0.4000 0.4157 0.4353 0.4510 0.4706 0.4863 0.5020 0.5059 0.5255 0.5451 0.5608 0.5765 0.5961 0.6118 0.6275 0.6353 0.6510 0.6706 0.6902 0.7020 0.7216 0.7373 0.7451 0.7608 0.7804 0.8000 0.8157 0.8275 0.8471 0.8510 0.8706 0.8902 0.9059 0.9255 0.9412 0.9569 0.9608 0.9804 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9922 0.9804 0.9765 0.9608 0.9608 0.9451 0.9451 0.9373 0.9255 0.9216 0.9059 0.9059 0.8902 0.8902 0.8824 0.8706 0.8667 0.8510 0.8510 0.8353 0.8353 0.8275 0.8157 0.8118 0.8000 0.7961 0.7804 0.7804 0.7608 0.7608 0.7569 0.7451 0.7412 0.7255 0.7255 0.7059 0.7059 0.7020 0.6902 0.6863 0.6706 0.6706 0.6510 0.6510 0.6353 0.6353 0.6314 0.6157 0.6157 0.6000 0.6000 0.5922 0.5804 0.5765 0.5608 0.5608 0.5451 0.5451 0.5373 0.5255 0.5059 0.5059 0.5059];
    case 'circular'
        red=[ 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9922 0.9765 0.9608 0.9451 0.9255 0.9059 0.9020 0.8863 0.8706 0.8510 0.8353 0.8157 0.8000 0.7922 0.7765 0.7608 0.7412 0.7255 0.7059 0.6902 0.6706 0.6510 0.6353 0.6275 0.6118 0.5961 0.5804 0.5608 0.5451 0.5255 0.5059 0.4902 0.4824 0.4667 0.4510 0.4353 0.4157 0.4000 0.3804 0.3765 0.3608 0.3451 0.3255 0.3059 0.3020 0.2902 0.2706 0.2510 0.2510 0.2353 0.2157 0.2118 0.2000 0.1804 0.1765 0.1608 0.1569 0.1451 0.1255 0.1255 0.1059 0.1059 0.0902 0.0902 0.0824 0.0706 0.0667 0.0510 0.0510 0.0471 0.0353 0.0353 0.0275 0.0157 0.0157 0.0157 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0118 0.0157 0.0157 0.0157 0.0275 0.0353 0.0353 0.0471 0.0510 0.0510 0.0667 0.0706 0.0824 0.0902 0.0902 0.1059 0.1059 0.1255 0.1255 0.1451 0.1569 0.1608 0.1765 0.1804 0.2000 0.2118 0.2157 0.2353 0.2510 0.2510 0.2706 0.2902 0.3020 0.3059 0.3255 0.3451 0.3608 0.3765 0.3804 0.4000 0.4157 0.4353 0.4510 0.4667 0.4824 0.4902 0.5059 0.5255 0.5451 0.5608 0.5804 0.5961 0.6118 0.6275 0.6353 0.6510 0.6706 0.6902 0.7059 0.7255 0.7412 0.7608 0.7765 0.7922 0.8000 0.8157 0.8353 0.8510 0.8706 0.8863 0.9020 0.9059 0.9255 0.9451 0.9608 0.9765 0.9922 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
        green=[ 0.3255 0.3059 0.3020 0.2902 0.2706 0.2510 0.2510 0.2353 0.2157 0.2118 0.2000 0.1804 0.1765 0.1608 0.1569 0.1451 0.1255 0.1255 0.1059 0.1059 0.0902 0.0902 0.0824 0.0706 0.0667 0.0510 0.0510 0.0471 0.0353 0.0353 0.0275 0.0157 0.0157 0.0157 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0118 0.0157 0.0157 0.0157 0.0275 0.0353 0.0353 0.0471 0.0510 0.0510 0.0667 0.0706 0.0824 0.0902 0.0902 0.1059 0.1059 0.1255 0.1255 0.1451 0.1569 0.1608 0.1765 0.1804 0.2000 0.2118 0.2157 0.2353 0.2510 0.2510 0.2706 0.2902 0.3020 0.3059 0.3255 0.3451 0.3608 0.3765 0.3804 0.4000 0.4157 0.4353 0.4510 0.4667 0.4824 0.4902 0.5059 0.5255 0.5451 0.5608 0.5804 0.5961 0.6118 0.6275 0.6353 0.6510 0.6706 0.6902 0.7059 0.7255 0.7412 0.7608 0.7765 0.7922 0.8000 0.8157 0.8353 0.8510 0.8706 0.8863 0.9020 0.9059 0.9255 0.9451 0.9608 0.9765 0.9922 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9922 0.9765 0.9608 0.9451 0.9255 0.9059 0.9020 0.8863 0.8706 0.8510 0.8353 0.8157 0.8000 0.7922 0.7765 0.7608 0.7412 0.7255 0.7059 0.6902 0.6706 0.6510 0.6353 0.6275 0.6118 0.5961 0.5804 0.5608 0.5451 0.5255 0.5059 0.4902 0.4824 0.4667 0.4510 0.4353 0.4157 0.4000 0.3804 0.3765 0.3608 0.3451 0.3255];
        blue=[ 0.3255 0.3451 0.3608 0.3765 0.3804 0.4000 0.4157 0.4353 0.4510 0.4667 0.4824 0.4902 0.5059 0.5255 0.5451 0.5608 0.5804 0.5961 0.6118 0.6275 0.6353 0.6510 0.6706 0.6902 0.7059 0.7255 0.7412 0.7608 0.7765 0.7922 0.8000 0.8157 0.8353 0.8510 0.8706 0.8863 0.9020 0.9059 0.9255 0.9451 0.9608 0.9765 0.9922 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9922 0.9765 0.9608 0.9451 0.9255 0.9059 0.9020 0.8863 0.8706 0.8510 0.8353 0.8157 0.8000 0.7922 0.7765 0.7608 0.7412 0.7255 0.7059 0.6902 0.6706 0.6510 0.6353 0.6275 0.6118 0.5961 0.5804 0.5608 0.5451 0.5255 0.5059 0.4902 0.4824 0.4667 0.4510 0.4353 0.4157 0.4000 0.3804 0.3765 0.3608 0.3451 0.3255 0.3059 0.3020 0.2902 0.2706 0.2510 0.2510 0.2353 0.2157 0.2118 0.2000 0.1804 0.1765 0.1608 0.1569 0.1451 0.1255 0.1255 0.1059 0.1059 0.0902 0.0902 0.0824 0.0706 0.0667 0.0510 0.0510 0.0471 0.0353 0.0353 0.0275 0.0157 0.0157 0.0157 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0118 0.0157 0.0157 0.0157 0.0275 0.0353 0.0353 0.0471 0.0510 0.0510 0.0667 0.0706 0.0824 0.0902 0.0902 0.1059 0.1059 0.1255 0.1255 0.1451 0.1569 0.1608 0.1765 0.1804 0.2000 0.2118 0.2157 0.2353 0.2510 0.2510 0.2706 0.2902 0.3020 0.3059 0.3255];
    case 'jfm base'
        red=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0118 0.0118 0.0118 0.0118 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0275 0.0275 0.0275 0.0314 0.0314 0.0353 0.0353 0.0353 0.0353 0.0353 0.0353 0.0353 0.0353 0.0471 0.0471 0.0510 0.0510 0.0510 0.0510 0.0510 0.0510 0.0510 0.0510 0.0667 0.0667 0.0706 0.0706 0.0706 0.0706 0.0706 0.0706 0.0824 0.0863 0.0863 0.0902 0.0902 0.0902 0.0902 0.0902 0.1020 0.1020 0.1059 0.1059 0.1059 0.1059 0.1059 0.1216 0.1216 0.1255 0.1255 0.1255 0.1255 0.1373 0.1412 0.1412 0.1451 0.1451 0.1451 0.1451 0.1569 0.1608 0.1608 0.1608 0.1608 0.1608 0.1765 0.1804 0.1804 0.1804 0.1922 0.1922 0.1961 0.2000 0.2000 0.2000 0.2118 0.2157 0.2157 0.2157 0.2275 0.2314 0.2353 0.2353 0.2353 0.2353 0.2471 0.2510 0.2510 0.2510 0.2510 0.2667 0.2706 0.2706 0.2824 0.2863 0.2902 0.2902 0.2902 0.3020 0.3059 0.3059 0.3059 0.3059 0.3216 0.3255 0.3255 0.3373 0.3412 0.3451 0.3451 0.3451 0.3569 0.3608 0.3608 0.3608 0.3765 0.3804 0.3804 0.3804 0.3961 0.4000 0.4000 0.4000 0.4118 0.4157 0.4157 0.4275 0.4314 0.4353 0.4353 0.4471 0.4510 0.4510 0.4510 0.4667 0.4706 0.4706 0.4706 0.4863 0.4902 0.4902 0.4902 0.5059 0.5059 0.5059 0.5059 0.5255 0.5255 0.5255 0.5373 0.5451 0.5451 0.5451 0.5569 0.5608 0.5608 0.5608 0.5765 0.5804 0.5804 0.5922 0.5961 0.6000 0.6000 0.6000 0.6118 0.6157 0.6157 0.6275 0.6314 0.6353 0.6353 0.6353 0.6471 0.6510 0.6510 0.6510 0.6667 0.6706 0.6706 0.6706 0.6706 0.6824 0.6863 0.6902 0.6902 0.6902 0.7020 0.7059 0.7059 0.7059 0.7059 0.7059 0.7059 0.7216 0.7255 0.7255 0.7255 0.7255 0.7255 0.7255 0.7373 0.7373 0.7412 0.7412 0.7412 0.7451 0.7451 0.7451 0.7451 0.7451 0.7451 0.7451];
        green=[ 0.0902 0.0902 0.0902 0.0902 0.0902 0.0902 0.0902 0.0902 0.0902 0.0902 0.0902 0.1020 0.1020 0.1020 0.1020 0.1020 0.1059 0.1059 0.1059 0.1059 0.1059 0.1059 0.1059 0.1059 0.1059 0.1059 0.1059 0.1216 0.1216 0.1216 0.1255 0.1255 0.1255 0.1255 0.1255 0.1255 0.1373 0.1412 0.1412 0.1451 0.1451 0.1451 0.1451 0.1451 0.1569 0.1569 0.1608 0.1608 0.1608 0.1608 0.1608 0.1765 0.1765 0.1804 0.1804 0.1804 0.1804 0.1922 0.1961 0.2000 0.2000 0.2000 0.2000 0.2118 0.2157 0.2157 0.2157 0.2275 0.2314 0.2314 0.2353 0.2353 0.2353 0.2471 0.2510 0.2510 0.2510 0.2510 0.2667 0.2706 0.2706 0.2706 0.2824 0.2863 0.2902 0.2902 0.2902 0.3020 0.3059 0.3059 0.3059 0.3059 0.3216 0.3255 0.3255 0.3255 0.3373 0.3412 0.3451 0.3451 0.3569 0.3608 0.3608 0.3608 0.3608 0.3765 0.3804 0.3804 0.3804 0.3961 0.4000 0.4000 0.4000 0.4118 0.4157 0.4157 0.4157 0.4314 0.4353 0.4353 0.4353 0.4471 0.4510 0.4510 0.4510 0.4667 0.4706 0.4706 0.4706 0.4824 0.4863 0.4902 0.4902 0.5020 0.5059 0.5059 0.5059 0.5059 0.5216 0.5255 0.5255 0.5373 0.5412 0.5451 0.5451 0.5451 0.5569 0.5608 0.5608 0.5608 0.5765 0.5804 0.5804 0.5804 0.5922 0.5961 0.6000 0.6000 0.6118 0.6157 0.6157 0.6157 0.6275 0.6314 0.6353 0.6353 0.6353 0.6471 0.6510 0.6510 0.6510 0.6510 0.6667 0.6706 0.6706 0.6706 0.6824 0.6863 0.6902 0.6902 0.6902 0.7020 0.7059 0.7059 0.7059 0.7059 0.7216 0.7255 0.7255 0.7255 0.7255 0.7373 0.7412 0.7451 0.7451 0.7451 0.7451 0.7569 0.7608 0.7608 0.7608 0.7608 0.7608 0.7765 0.7765 0.7804 0.7804 0.7804 0.7804 0.7922 0.7961 0.7961 0.8000 0.8000 0.8000 0.8000 0.8000 0.8118 0.8118 0.8157 0.8157 0.8157 0.8157 0.8157 0.8275 0.8275 0.8314 0.8314 0.8314 0.8353 0.8353 0.8353 0.8353 0.8353 0.8353 0.8353 0.8353 0.8471 0.8471 0.8471 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510];
        blue=[ 0.2706 0.2706 0.2706 0.2863 0.2902 0.2902 0.3020 0.3059 0.3059 0.3059 0.3216 0.3255 0.3255 0.3373 0.3412 0.3451 0.3451 0.3569 0.3608 0.3608 0.3608 0.3765 0.3804 0.3804 0.3922 0.3961 0.4000 0.4000 0.4118 0.4157 0.4157 0.4157 0.4314 0.4353 0.4353 0.4353 0.4471 0.4510 0.4510 0.4510 0.4667 0.4706 0.4706 0.4706 0.4824 0.4902 0.4902 0.4902 0.5020 0.5059 0.5059 0.5059 0.5059 0.5216 0.5255 0.5255 0.5373 0.5412 0.5451 0.5451 0.5451 0.5569 0.5608 0.5608 0.5608 0.5608 0.5765 0.5804 0.5804 0.5804 0.5922 0.5961 0.6000 0.6000 0.6000 0.6118 0.6157 0.6157 0.6157 0.6275 0.6314 0.6353 0.6353 0.6353 0.6353 0.6471 0.6510 0.6510 0.6510 0.6510 0.6667 0.6706 0.6706 0.6706 0.6706 0.6824 0.6863 0.6902 0.6902 0.6902 0.6902 0.7020 0.7059 0.7059 0.7059 0.7059 0.7059 0.7216 0.7255 0.7255 0.7255 0.7255 0.7373 0.7412 0.7412 0.7451 0.7451 0.7451 0.7451 0.7569 0.7608 0.7608 0.7608 0.7608 0.7608 0.7608 0.7765 0.7804 0.7804 0.7804 0.7804 0.7804 0.7922 0.7961 0.7961 0.8000 0.8000 0.8000 0.8000 0.8118 0.8118 0.8157 0.8157 0.8157 0.8157 0.8157 0.8275 0.8275 0.8314 0.8314 0.8353 0.8353 0.8353 0.8353 0.8353 0.8471 0.8471 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8510 0.8667 0.8667 0.8706 0.8706 0.8706 0.8706 0.8706 0.8706 0.8824 0.8824 0.8824 0.8863 0.8863 0.8902 0.8902 0.8902 0.8902 0.8902 0.8902 0.8902 0.8902 0.9020 0.9020 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9216 0.9216 0.9216 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9373 0.9373 0.9373 0.9373 0.9373 0.9412 0.9412 0.9412 0.9412 0.9412 0.9412 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451];
    case 'jfm extended'
        red=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0118 0.0118 0.0118 0.0118 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0275 0.0275 0.0275 0.0314 0.0314 0.0353 0.0353 0.0353 0.0353 0.0353 0.0353 0.0471 0.0471 0.0510 0.0510 0.0510 0.0510 0.0510 0.0510 0.0510 0.0510 0.0667 0.0667 0.0706 0.0706 0.0706 0.0706 0.0824 0.0824 0.0863 0.0902 0.0902 0.0902 0.0902 0.0902 0.1020 0.1059 0.1059 0.1059 0.1059 0.1059 0.1216 0.1216 0.1255 0.1255 0.1255 0.1373 0.1412 0.1412 0.1451 0.1451 0.1451 0.1569 0.1608 0.1608 0.1608 0.1608 0.1765 0.1804 0.1804 0.1804 0.1804 0.1922 0.1961 0.2000 0.2000 0.2118 0.2157 0.2157 0.2157 0.2275 0.2314 0.2353 0.2353 0.2353 0.2471 0.2510 0.2510 0.2510 0.2667 0.2706 0.2706 0.2706 0.2824 0.2902 0.2902 0.2902 0.3020 0.3059 0.3059 0.3059 0.3216 0.3255 0.3255 0.3373 0.3412 0.3451 0.3451 0.3569 0.3608 0.3608 0.3608 0.3804 0.3804 0.3922 0.3961 0.4000 0.4000 0.4118 0.4157 0.4157 0.4275 0.4353 0.4353 0.4471 0.4510 0.4510 0.4510 0.4706 0.4706 0.4824 0.4863 0.4902 0.4902 0.5059 0.5059 0.5059 0.5216 0.5255 0.5255 0.5412 0.5451 0.5451 0.5569 0.5608 0.5608 0.5765 0.5804 0.5804 0.5961 0.6000 0.6000 0.6118 0.6157 0.6275 0.6314 0.6353 0.6353 0.6510 0.6510 0.6510 0.6706 0.6706 0.6824 0.6863 0.6902 0.7020 0.7059 0.7059 0.7059 0.7255 0.7255 0.7373 0.7412 0.7451 0.7451 0.7608 0.7608 0.7608 0.7765 0.7804 0.7804 0.7961 0.8000 0.8000 0.8118 0.8157 0.8157 0.8275 0.8353 0.8353 0.8353 0.8471 0.8510 0.8510 0.8510 0.8706 0.8706 0.8706 0.8824 0.8863 0.8902 0.8902 0.8902 0.9020 0.9059 0.9059 0.9059 0.9059 0.9216 0.9255 0.9255 0.9255 0.9255 0.9373 0.9373 0.9412 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451];
        green=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0118 0.0118 0.0118 0.0157 0.0157 0.0157 0.0157 0.0157 0.0157 0.0275 0.0275 0.0314 0.0314 0.0353 0.0353 0.0353 0.0353 0.0353 0.0471 0.0510 0.0510 0.0510 0.0510 0.0510 0.0667 0.0667 0.0706 0.0706 0.0706 0.0824 0.0824 0.0863 0.0902 0.0902 0.0902 0.1020 0.1059 0.1059 0.1059 0.1059 0.1216 0.1255 0.1255 0.1255 0.1373 0.1412 0.1451 0.1451 0.1451 0.1569 0.1608 0.1608 0.1608 0.1765 0.1804 0.1804 0.1804 0.1961 0.2000 0.2000 0.2000 0.2118 0.2157 0.2157 0.2275 0.2314 0.2353 0.2353 0.2471 0.2510 0.2510 0.2510 0.2667 0.2706 0.2706 0.2824 0.2902 0.2902 0.2902 0.3059 0.3059 0.3059 0.3216 0.3255 0.3255 0.3373 0.3412 0.3451 0.3451 0.3608 0.3608 0.3608 0.3765 0.3804 0.3804 0.3922 0.3961 0.4000 0.4000 0.4157 0.4157 0.4275 0.4314 0.4353 0.4353 0.4471 0.4510 0.4510 0.4667 0.4706 0.4706 0.4824 0.4863 0.4902 0.4902 0.5059 0.5059 0.5059 0.5216 0.5255 0.5255 0.5373 0.5451 0.5451 0.5569 0.5608 0.5608 0.5608 0.5765 0.5804 0.5804 0.5961 0.6000 0.6000 0.6118 0.6157 0.6157 0.6275 0.6314 0.6353 0.6353 0.6471 0.6510 0.6510 0.6667 0.6706 0.6706 0.6824 0.6863 0.6902 0.6902 0.7020 0.7059 0.7059 0.7059 0.7216 0.7255 0.7255 0.7373 0.7412 0.7451 0.7451 0.7569 0.7608 0.7608 0.7608 0.7608 0.7804 0.7804 0.7804 0.7922 0.7961 0.8000 0.8000 0.8000 0.8157 0.8157 0.8157 0.8275 0.8314 0.8353 0.8353 0.8353 0.8471 0.8510 0.8510 0.8510 0.8510 0.8510 0.8667 0.8706 0.8706 0.8706 0.8824 0.8863 0.8863 0.8902 0.8902 0.8902 0.9020 0.9020 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9216 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9373 0.9373 0.9412 0.9412 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9569 0.9569 0.9569 0.9569 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608];
        blue=[ 0.0000 0.0000 0.0118 0.0157 0.0275 0.0353 0.0353 0.0510 0.0510 0.0667 0.0706 0.0824 0.0902 0.0902 0.1059 0.1059 0.1059 0.1255 0.1255 0.1412 0.1451 0.1569 0.1608 0.1608 0.1765 0.1804 0.1922 0.2000 0.2000 0.2118 0.2157 0.2275 0.2314 0.2353 0.2471 0.2510 0.2510 0.2667 0.2706 0.2824 0.2863 0.2902 0.3020 0.3059 0.3059 0.3216 0.3255 0.3255 0.3373 0.3451 0.3451 0.3569 0.3608 0.3608 0.3765 0.3804 0.3804 0.3961 0.4000 0.4000 0.4118 0.4157 0.4157 0.4314 0.4353 0.4353 0.4471 0.4510 0.4510 0.4667 0.4706 0.4706 0.4824 0.4902 0.4902 0.4902 0.5059 0.5059 0.5059 0.5216 0.5255 0.5255 0.5373 0.5412 0.5451 0.5451 0.5569 0.5608 0.5608 0.5608 0.5804 0.5804 0.5804 0.5961 0.6000 0.6000 0.6000 0.6157 0.6157 0.6157 0.6275 0.6314 0.6353 0.6353 0.6471 0.6510 0.6510 0.6510 0.6667 0.6706 0.6706 0.6706 0.6824 0.6863 0.6902 0.6902 0.7020 0.7059 0.7059 0.7059 0.7059 0.7216 0.7255 0.7255 0.7255 0.7373 0.7412 0.7451 0.7451 0.7451 0.7569 0.7608 0.7608 0.7608 0.7608 0.7765 0.7804 0.7804 0.7804 0.7922 0.7961 0.8000 0.8000 0.8000 0.8000 0.8118 0.8157 0.8157 0.8157 0.8157 0.8275 0.8314 0.8353 0.8353 0.8353 0.8353 0.8471 0.8510 0.8510 0.8510 0.8510 0.8510 0.8667 0.8667 0.8706 0.8706 0.8706 0.8706 0.8824 0.8824 0.8863 0.8902 0.8902 0.8902 0.8902 0.8902 0.9020 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9059 0.9216 0.9216 0.9255 0.9255 0.9255 0.9255 0.9255 0.9255 0.9373 0.9373 0.9412 0.9412 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9451 0.9569 0.9569 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9608 0.9765 0.9765 0.9765 0.9765 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9804 0.9922 0.9922 0.9922 0.9922 0.9922 0.9922 0.9922 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 0.9961 1.0000];
    case 'bipolar'
        red=[ 0.4980 0.4902 0.4824 0.4745 0.4667 0.4588 0.4510 0.4431 0.4353 0.4275 0.4196 0.4118 0.4039 0.3961 0.3882 0.3804 0.3725 0.3647 0.3569 0.3490 0.3412 0.3333 0.3255 0.3176 0.3098 0.3020 0.2941 0.2863 0.2784 0.2706 0.2627 0.2549 0.2471 0.2392 0.2314 0.2235 0.2157 0.2078 0.2000 0.1922 0.1843 0.1765 0.1686 0.1608 0.1529 0.1451 0.1373 0.1294 0.1216 0.1137 0.1059 0.0980 0.0902 0.0824 0.0745 0.0667 0.0588 0.0510 0.0431 0.0353 0.0275 0.0196 0.0118 0.0039 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
        green=[ 1.0000 0.9843 0.9686 0.9529 0.9373 0.9216 0.9059 0.8902 0.8745 0.8588 0.8431 0.8275 0.8118 0.7961 0.7804 0.7647 0.7490 0.7333 0.7176 0.7020 0.6863 0.6706 0.6549 0.6392 0.6235 0.6078 0.5922 0.5765 0.5608 0.5451 0.5294 0.5137 0.4980 0.4824 0.4667 0.4510 0.4353 0.4196 0.4039 0.3882 0.3725 0.3569 0.3412 0.3255 0.3098 0.2941 0.2784 0.2627 0.2471 0.2314 0.2157 0.2000 0.1843 0.1686 0.1529 0.1373 0.1216 0.1059 0.0902 0.0745 0.0588 0.0431 0.0275 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882];
        blue=[ 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9843 0.9686 0.9529 0.9373 0.9216 0.9059 0.8902 0.8745 0.8588 0.8431 0.8275 0.8118 0.7961 0.7804 0.7647 0.7490 0.7333 0.7176 0.7020 0.6863 0.6706 0.6549 0.6392 0.6235 0.6078 0.5922 0.5765 0.5608 0.5451 0.5294 0.5137 0.4980 0.4824 0.4667 0.4510 0.4353 0.4196 0.4039 0.3882 0.3725 0.3569 0.3412 0.3255 0.3098 0.2941 0.2784 0.2627 0.2471 0.2314 0.2157 0.2000 0.1843 0.1686 0.1529 0.1373 0.1216 0.1059 0.0902 0.0745 0.0588 0.0431 0.0275 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0078 0.0157 0.0235 0.0314 0.0392 0.0471 0.0549 0.0627 0.0706 0.0784 0.0863 0.0941 0.1020 0.1098 0.1176 0.1255 0.1333 0.1412 0.1490 0.1569 0.1647 0.1725 0.1804 0.1882 0.1961 0.2039 0.2118 0.2196 0.2275 0.2353 0.2431 0.2510 0.2588 0.2667 0.2745 0.2824 0.2902 0.2980 0.3059 0.3137 0.3216 0.3294 0.3373 0.3451 0.3529 0.3608 0.3686 0.3765 0.3843 0.3922 0.4000 0.4078 0.4157 0.4235 0.4314 0.4392 0.4471 0.4549 0.4627 0.4706 0.4784 0.4863 0.4941];
    case 'schlieren'
        red=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0275 0.0392 0.0510 0.0627 0.0745 0.0863 0.0980 0.1098 0.1216 0.1333 0.1451 0.1569 0.1686 0.1804 0.1922 0.2039 0.2157 0.2275 0.2392 0.2510 0.2627 0.2745 0.2863 0.2980 0.3098 0.3216 0.3333 0.3451 0.3569 0.3686 0.3804 0.3922 0.4039 0.4157 0.4275 0.4392 0.4510 0.4627 0.4745 0.4863 0.4980 0.5098 0.5216 0.5333 0.5451 0.5569 0.5686 0.5804 0.5922 0.6039 0.6157 0.6275 0.6392 0.6510 0.6627 0.6745 0.6863 0.6980 0.7098 0.7216 0.7333 0.7451 0.7569 0.7686 0.7804 0.7922 0.8039 0.8157 0.8275 0.8392 0.8510 0.8627 0.8745 0.8863 0.8980 0.9098 0.9216 0.9333 0.9451 0.9569 0.9686 0.9804 0.9922 1.0000];
        green=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0275 0.0392 0.0510 0.0627 0.0745 0.0863 0.0980 0.1098 0.1216 0.1333 0.1451 0.1569 0.1686 0.1804 0.1922 0.2039 0.2157 0.2275 0.2392 0.2510 0.2627 0.2745 0.2863 0.2980 0.3098 0.3216 0.3333 0.3451 0.3569 0.3686 0.3804 0.3922 0.4039 0.4157 0.4275 0.4392 0.4510 0.4627 0.4745 0.4863 0.4980 0.5098 0.5216 0.5333 0.5451 0.5569 0.5686 0.5804 0.5922 0.6039 0.6157 0.6275 0.6392 0.6510 0.6627 0.6745 0.6863 0.6980 0.7098 0.7216 0.7333 0.7451 0.7569 0.7686 0.7804 0.7922 0.8039 0.8157 0.8275 0.8392 0.8510 0.8627 0.8745 0.8863 0.8980 0.9098 0.9216 0.9333 0.9451 0.9569 0.9686 0.9804 0.9922 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
        blue=[ 0.0000 0.0118 0.0235 0.0353 0.0471 0.0588 0.0706 0.0824 0.0941 0.1059 0.1176 0.1294 0.1412 0.1529 0.1647 0.1765 0.1882 0.2000 0.2118 0.2235 0.2353 0.2471 0.2588 0.2706 0.2824 0.2941 0.3059 0.3176 0.3294 0.3412 0.3529 0.3647 0.3765 0.3882 0.4000 0.4118 0.4235 0.4353 0.4471 0.4588 0.4706 0.4824 0.4941 0.5059 0.5176 0.5294 0.5412 0.5529 0.5647 0.5765 0.5882 0.6000 0.6118 0.6235 0.6353 0.6471 0.6588 0.6706 0.6824 0.6941 0.7059 0.7176 0.7294 0.7412 0.7529 0.7647 0.7765 0.7882 0.8000 0.8118 0.8235 0.8353 0.8471 0.8588 0.8706 0.8824 0.8941 0.9059 0.9176 0.9294 0.9412 0.9529 0.9647 0.9765 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
    case 'single cycle - double brightness'
        red=[ 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9843 0.9686 0.9529 0.9373 0.9216 0.9059 0.8902 0.8745 0.8588 0.8431 0.8275 0.8118 0.7961 0.7804 0.7647 0.7490 0.7333 0.7176 0.7020 0.6863 0.6706 0.6549 0.6392 0.6235 0.6078 0.5922 0.5765 0.5608 0.5451 0.5294 0.5137 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
        green=[ 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9843 0.9686 0.9529 0.9373 0.9216 0.9059 0.8902 0.8745 0.8588 0.8431 0.8275 0.8118 0.7961 0.7804 0.7647 0.7490 0.7333 0.7176 0.7020 0.6863 0.6706 0.6549 0.6392 0.6235 0.6078 0.5922 0.5765 0.5608 0.5451 0.5294 0.5137 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5098 0.5176 0.5255 0.5333 0.5412 0.5490 0.5569 0.5647 0.5725 0.5804 0.5882 0.5961 0.6039 0.6118 0.6196 0.6275 0.6353 0.6431 0.6510 0.6588 0.6667 0.6745 0.6824 0.6902 0.6980 0.7059 0.7137 0.7216 0.7294 0.7373 0.7451 0.7529 0.7608 0.7686 0.7765 0.7843 0.7922 0.8000 0.8078 0.8157 0.8235 0.8314 0.8392 0.8471 0.8549 0.8627 0.8706 0.8784 0.8863 0.8941 0.9020 0.9098 0.9176 0.9255 0.9333 0.9412 0.9490 0.9569 0.9647 0.9725 0.9804 0.9882 0.9961];
        blue=[ 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9922 0.9843 0.9765 0.9686 0.9608 0.9529 0.9451 0.9373 0.9294 0.9216 0.9137 0.9059 0.8980 0.8902 0.8824 0.8745 0.8667 0.8588 0.8510 0.8431 0.8353 0.8275 0.8196 0.8118 0.8039 0.7961 0.7882 0.7804 0.7725 0.7647 0.7569 0.7529 0.7608 0.7686 0.7765 0.7843 0.7922 0.8000 0.8078 0.8157 0.8235 0.8314 0.8392 0.8471 0.8549 0.8627 0.8706 0.8784 0.8863 0.8941 0.9020 0.9098 0.9176 0.9255 0.9333 0.9412 0.9490 0.9569 0.9647 0.9725 0.9804 0.9882 0.9961];
    case 'single cycle - aperture'
        red=[ 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4824 0.4667 0.4510 0.4353 0.4196 0.4039 0.3882 0.3725 0.3569 0.3412 0.3255 0.3098 0.2941 0.2784 0.2627 0.2471 0.2314 0.2157 0.2000 0.1843 0.1686 0.1529 0.1373 0.1216 0.1059 0.0902 0.0745 0.0588 0.0431 0.0275 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 1.0000];
        green=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4824 0.4667 0.4510 0.4353 0.4196 0.4039 0.3882 0.3725 0.3569 0.3412 0.3255 0.3098 0.2941 0.2784 0.2627 0.2471 0.2314 0.2157 0.2000 0.1843 0.1686 0.1529 0.1373 0.1216 0.1059 0.0902 0.0745 0.0588 0.0431 0.0275 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0078 0.0157 0.0235 0.0314 0.0392 0.0471 0.0549 0.0627 0.0706 0.0784 0.0863 0.0941 0.1020 0.1098 0.1176 0.1255 0.1333 0.1412 0.1490 0.1569 0.1647 0.1725 0.1804 0.1882 0.1961 0.2039 0.2118 0.2196 0.2275 0.2353 0.2431 0.2510 0.2588 0.2667 0.2745 0.2824 0.2902 0.2980 0.3059 0.3137 0.3216 0.3294 0.3373 0.3451 0.3529 0.3608 0.3686 0.3765 0.3843 0.3922 0.4000 0.4078 0.4157 0.4235 0.4314 0.4392 0.4471 0.4549 0.4627 0.4706 0.4784 0.4863 1.0000];
        blue=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4902 0.4824 0.4745 0.4667 0.4588 0.4510 0.4431 0.4353 0.4275 0.4196 0.4118 0.4039 0.3961 0.3882 0.3804 0.3725 0.3647 0.3569 0.3490 0.3412 0.3333 0.3255 0.3176 0.3098 0.3020 0.2941 0.2863 0.2784 0.2706 0.2627 0.2549 0.2510 0.2588 0.2667 0.2745 0.2824 0.2902 0.2980 0.3059 0.3137 0.3216 0.3294 0.3373 0.3451 0.3529 0.3608 0.3686 0.3765 0.3843 0.3922 0.4000 0.4078 0.4157 0.4235 0.4314 0.4392 0.4471 0.4549 0.4627 0.4706 0.4784 0.4863 1.0000];
    case 'single cycle - half brightness'
        red=[ 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4824 0.4667 0.4510 0.4353 0.4196 0.4039 0.3882 0.3725 0.3569 0.3412 0.3255 0.3098 0.2941 0.2784 0.2627 0.2471 0.2314 0.2157 0.2000 0.1843 0.1686 0.1529 0.1373 0.1216 0.1059 0.0902 0.0745 0.0588 0.0431 0.0275 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980];
        green=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4824 0.4667 0.4510 0.4353 0.4196 0.4039 0.3882 0.3725 0.3569 0.3412 0.3255 0.3098 0.2941 0.2784 0.2627 0.2471 0.2314 0.2157 0.2000 0.1843 0.1686 0.1529 0.1373 0.1216 0.1059 0.0902 0.0745 0.0588 0.0431 0.0275 0.0118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0078 0.0157 0.0235 0.0314 0.0392 0.0471 0.0549 0.0627 0.0706 0.0784 0.0863 0.0941 0.1020 0.1098 0.1176 0.1255 0.1333 0.1412 0.1490 0.1569 0.1647 0.1725 0.1804 0.1882 0.1961 0.2039 0.2118 0.2196 0.2275 0.2353 0.2431 0.2510 0.2588 0.2667 0.2745 0.2824 0.2902 0.2980 0.3059 0.3137 0.3216 0.3294 0.3373 0.3451 0.3529 0.3608 0.3686 0.3765 0.3843 0.3922 0.4000 0.4078 0.4157 0.4235 0.4314 0.4392 0.4471 0.4549 0.4627 0.4706 0.4784 0.4863 0.4941];
        blue=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4980 0.4902 0.4824 0.4745 0.4667 0.4588 0.4510 0.4431 0.4353 0.4275 0.4196 0.4118 0.4039 0.3961 0.3882 0.3804 0.3725 0.3647 0.3569 0.3490 0.3412 0.3333 0.3255 0.3176 0.3098 0.3020 0.2941 0.2863 0.2784 0.2706 0.2627 0.2549 0.2510 0.2588 0.2667 0.2745 0.2824 0.2902 0.2980 0.3059 0.3137 0.3216 0.3294 0.3373 0.3451 0.3529 0.3608 0.3686 0.3765 0.3843 0.3922 0.4000 0.4078 0.4157 0.4235 0.4314 0.4392 0.4471 0.4549 0.4627 0.4706 0.4784 0.4863 0.4941];
    case '(default)'
        red=[ 0.0000 0.0314 0.0510 0.0902 0.1255 0.1569 0.1804 0.2157 0.2510 0.2824 0.3059 0.3451 0.3765 0.4000 0.4353 0.4706 0.5020 0.5255 0.5608 0.5961 0.6275 0.6510 0.6902 0.7216 0.7451 0.7804 0.8157 0.8471 0.8706 0.9059 0.9412 0.9608 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9608 0.9373 0.9059 0.8706 0.8353 0.8118 0.7804 0.7451 0.7059 0.6863 0.6510 0.6157 0.5922 0.5608 0.5255 0.4902 0.4667 0.4353 0.4000 0.3608 0.3412 0.3059 0.2706 0.2471 0.2157 0.1804 0.1451 0.1216 0.0902 0.0510 0.0275 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0314 0.0510 0.0902 0.1255 0.1569 0.1804 0.2157 0.2510 0.2824 0.3059 0.3451 0.3765 0.4000 0.4353 0.4706 0.5020 0.5255 0.5608 0.5961 0.6275 0.6510 0.6902 0.7216 0.7451 0.7804 0.8157 0.8471 0.8706 0.9059 0.9412 0.9608 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
        green=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0314 0.0510 0.0902 0.1255 0.1569 0.1804 0.2157 0.2510 0.2824 0.3059 0.3451 0.3765 0.4000 0.4353 0.4706 0.5020 0.5255 0.5608 0.5961 0.6275 0.6510 0.6902 0.7216 0.7451 0.7804 0.8157 0.8471 0.8706 0.9059 0.9412 0.9608 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9608 0.9373 0.9059 0.8706 0.8353 0.8118 0.7804 0.7451 0.7059 0.6863 0.6510 0.6157 0.5922 0.5608 0.5255 0.4902 0.4667 0.4353 0.4000 0.3608 0.3412 0.3059 0.2706 0.2471 0.2157 0.1804 0.1451 0.1216 0.0902 0.0510 0.0275 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0510 0.0706 0.0902 0.1059 0.1255 0.1412 0.1569 0.1608 0.1804 0.2000 0.2157 0.2353 0.2510 0.2667 0.2824 0.2902 0.3059 0.3255 0.3451 0.3608 0.3765 0.3922 0.4000 0.4157 0.4353 0.4510 0.4706 0.4863 0.5020 0.5059 0.5255 0.5451 0.5608 0.5804 0.5961 0.6118 0.6275 0.6353 0.6510 0.6706 0.6902 0.7059 0.7216 0.7373 0.7451 0.7608 0.7804 0.8000 0.8157 0.8314 0.8471 0.8510 0.8706 0.8902 0.9059 0.9255 0.9412 0.9569 0.9608 0.9804];
        blue=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0314 0.0510 0.0902 0.1255 0.1569 0.1804 0.2157 0.2510 0.2824 0.3059 0.3451 0.3765 0.4000 0.4353 0.4706 0.5020 0.5255 0.5608 0.5961 0.6275 0.6510 0.6902 0.7216 0.7451 0.7804 0.8157 0.8471 0.8706 0.9059 0.9412 0.9608 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9804 0.9608 0.9451 0.9373 0.9216 0.9059 0.8902 0.8706 0.8510 0.8353 0.8275 0.8118 0.7961 0.7804 0.7608 0.7451 0.7255 0.7059 0.7020 0.6863 0.6706 0.6510 0.6353 0.6157 0.6000 0.5922 0.5765 0.5608 0.5451 0.5255 0.5059 0.5020 0.5059 0.5255 0.5451 0.5608 0.5804 0.5961 0.6118 0.6275 0.6353 0.6510 0.6706 0.6902 0.7059 0.7216 0.7373 0.7451 0.7608 0.7804 0.8000 0.8157 0.8314 0.8471 0.8510 0.8706 0.8902 0.9059 0.9255 0.9412 0.9569 0.9608 0.9804];
    case 'single cycle'
        red=[ 0.0000 0.0314 0.0627 0.0941 0.1255 0.1569 0.1882 0.2196 0.2510 0.2824 0.3137 0.3451 0.3765 0.4078 0.4392 0.4706 0.5020 0.5333 0.5647 0.5961 0.6275 0.6588 0.6902 0.7216 0.7529 0.7843 0.8157 0.8471 0.8784 0.9098 0.9412 0.9725 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9686 0.9373 0.9059 0.8745 0.8431 0.8118 0.7804 0.7490 0.7176 0.6863 0.6549 0.6235 0.5922 0.5608 0.5294 0.4980 0.4667 0.4353 0.4039 0.3725 0.3412 0.3098 0.2784 0.2471 0.2157 0.1843 0.1529 0.1216 0.0902 0.0588 0.0275 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0314 0.0627 0.0941 0.1255 0.1569 0.1882 0.2196 0.2510 0.2824 0.3137 0.3451 0.3765 0.4078 0.4392 0.4706 0.5020 0.5333 0.5647 0.5961 0.6275 0.6588 0.6902 0.7216 0.7529 0.7843 0.8157 0.8471 0.8784 0.9098 0.9412 0.9725 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
        green=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0314 0.0627 0.0941 0.1255 0.1569 0.1882 0.2196 0.2510 0.2824 0.3137 0.3451 0.3765 0.4078 0.4392 0.4706 0.5020 0.5333 0.5647 0.5961 0.6275 0.6588 0.6902 0.7216 0.7529 0.7843 0.8157 0.8471 0.8784 0.9098 0.9412 0.9725 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9686 0.9373 0.9059 0.8745 0.8431 0.8118 0.7804 0.7490 0.7176 0.6863 0.6549 0.6235 0.5922 0.5608 0.5294 0.4980 0.4667 0.4353 0.4039 0.3725 0.3412 0.3098 0.2784 0.2471 0.2157 0.1843 0.1529 0.1216 0.0902 0.0588 0.0275 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0157 0.0314 0.0471 0.0627 0.0784 0.0941 0.1098 0.1255 0.1412 0.1569 0.1725 0.1882 0.2039 0.2196 0.2353 0.2510 0.2667 0.2824 0.2980 0.3137 0.3294 0.3451 0.3608 0.3765 0.3922 0.4078 0.4235 0.4392 0.4549 0.4706 0.4863 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882];
        blue=[ 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0314 0.0627 0.0941 0.1255 0.1569 0.1882 0.2196 0.2510 0.2824 0.3137 0.3451 0.3765 0.4078 0.4392 0.4706 0.5020 0.5333 0.5647 0.5961 0.6275 0.6588 0.6902 0.7216 0.7529 0.7843 0.8157 0.8471 0.8784 0.9098 0.9412 0.9725 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9843 0.9686 0.9529 0.9373 0.9216 0.9059 0.8902 0.8745 0.8588 0.8431 0.8275 0.8118 0.7961 0.7804 0.7647 0.7490 0.7333 0.7176 0.7020 0.6863 0.6706 0.6549 0.6392 0.6235 0.6078 0.5922 0.5765 0.5608 0.5451 0.5294 0.5137 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882];
    case 'negative'
        red=[ 1.0000 0.9961 0.9922 0.9882 0.9843 0.9804 0.9765 0.9725 0.9686 0.9647 0.9608 0.9569 0.9529 0.9490 0.9451 0.9412 0.9373 0.9333 0.9294 0.9255 0.9216 0.9176 0.9137 0.9098 0.9059 0.9020 0.8980 0.8941 0.8902 0.8863 0.8824 0.8784 0.8745 0.8706 0.8667 0.8627 0.8588 0.8549 0.8510 0.8471 0.8431 0.8392 0.8353 0.8314 0.8275 0.8235 0.8196 0.8157 0.8118 0.8078 0.8039 0.8000 0.7961 0.7922 0.7882 0.7843 0.7804 0.7765 0.7725 0.7686 0.7647 0.7608 0.7569 0.7529 0.7490 0.7451 0.7412 0.7373 0.7333 0.7294 0.7255 0.7216 0.7176 0.7137 0.7098 0.7059 0.7020 0.6980 0.6941 0.6902 0.6863 0.6824 0.6784 0.6745 0.6706 0.6667 0.6627 0.6588 0.6549 0.6510 0.6471 0.6431 0.6392 0.6353 0.6314 0.6275 0.6235 0.6196 0.6157 0.6118 0.6078 0.6039 0.6000 0.5961 0.5922 0.5882 0.5843 0.5804 0.5765 0.5725 0.5686 0.5647 0.5608 0.5569 0.5529 0.5490 0.5451 0.5412 0.5373 0.5333 0.5294 0.5255 0.5216 0.5176 0.5137 0.5098 0.5059 0.5020 0.4980 0.4941 0.4902 0.4863 0.4824 0.4784 0.4745 0.4706 0.4667 0.4627 0.4588 0.4549 0.4510 0.4471 0.4431 0.4392 0.4353 0.4314 0.4275 0.4235 0.4196 0.4157 0.4118 0.4078 0.4039 0.4000 0.3961 0.3922 0.3882 0.3843 0.3804 0.3765 0.3725 0.3686 0.3647 0.3608 0.3569 0.3529 0.3490 0.3451 0.3412 0.3373 0.3333 0.3294 0.3255 0.3216 0.3176 0.3137 0.3098 0.3059 0.3020 0.2980 0.2941 0.2902 0.2863 0.2824 0.2784 0.2745 0.2706 0.2667 0.2627 0.2588 0.2549 0.2510 0.2471 0.2431 0.2392 0.2353 0.2314 0.2275 0.2235 0.2196 0.2157 0.2118 0.2078 0.2039 0.2000 0.1961 0.1922 0.1882 0.1843 0.1804 0.1765 0.1725 0.1686 0.1647 0.1608 0.1569 0.1529 0.1490 0.1451 0.1412 0.1373 0.1333 0.1294 0.1255 0.1216 0.1176 0.1137 0.1098 0.1059 0.1020 0.0980 0.0941 0.0902 0.0863 0.0824 0.0784 0.0745 0.0706 0.0667 0.0627 0.0588 0.0549 0.0510 0.0471 0.0431 0.0392 0.0353 0.0314 0.0275 0.0235 0.0196 0.0157 0.0118 0.0078 0.0039 0.0000];
        green=[ 1.0000 0.9961 0.9922 0.9882 0.9843 0.9804 0.9765 0.9725 0.9686 0.9647 0.9608 0.9569 0.9529 0.9490 0.9451 0.9412 0.9373 0.9333 0.9294 0.9255 0.9216 0.9176 0.9137 0.9098 0.9059 0.9020 0.8980 0.8941 0.8902 0.8863 0.8824 0.8784 0.8745 0.8706 0.8667 0.8627 0.8588 0.8549 0.8510 0.8471 0.8431 0.8392 0.8353 0.8314 0.8275 0.8235 0.8196 0.8157 0.8118 0.8078 0.8039 0.8000 0.7961 0.7922 0.7882 0.7843 0.7804 0.7765 0.7725 0.7686 0.7647 0.7608 0.7569 0.7529 0.7490 0.7451 0.7412 0.7373 0.7333 0.7294 0.7255 0.7216 0.7176 0.7137 0.7098 0.7059 0.7020 0.6980 0.6941 0.6902 0.6863 0.6824 0.6784 0.6745 0.6706 0.6667 0.6627 0.6588 0.6549 0.6510 0.6471 0.6431 0.6392 0.6353 0.6314 0.6275 0.6235 0.6196 0.6157 0.6118 0.6078 0.6039 0.6000 0.5961 0.5922 0.5882 0.5843 0.5804 0.5765 0.5725 0.5686 0.5647 0.5608 0.5569 0.5529 0.5490 0.5451 0.5412 0.5373 0.5333 0.5294 0.5255 0.5216 0.5176 0.5137 0.5098 0.5059 0.5020 0.4980 0.4941 0.4902 0.4863 0.4824 0.4784 0.4745 0.4706 0.4667 0.4627 0.4588 0.4549 0.4510 0.4471 0.4431 0.4392 0.4353 0.4314 0.4275 0.4235 0.4196 0.4157 0.4118 0.4078 0.4039 0.4000 0.3961 0.3922 0.3882 0.3843 0.3804 0.3765 0.3725 0.3686 0.3647 0.3608 0.3569 0.3529 0.3490 0.3451 0.3412 0.3373 0.3333 0.3294 0.3255 0.3216 0.3176 0.3137 0.3098 0.3059 0.3020 0.2980 0.2941 0.2902 0.2863 0.2824 0.2784 0.2745 0.2706 0.2667 0.2627 0.2588 0.2549 0.2510 0.2471 0.2431 0.2392 0.2353 0.2314 0.2275 0.2235 0.2196 0.2157 0.2118 0.2078 0.2039 0.2000 0.1961 0.1922 0.1882 0.1843 0.1804 0.1765 0.1725 0.1686 0.1647 0.1608 0.1569 0.1529 0.1490 0.1451 0.1412 0.1373 0.1333 0.1294 0.1255 0.1216 0.1176 0.1137 0.1098 0.1059 0.1020 0.0980 0.0941 0.0902 0.0863 0.0824 0.0784 0.0745 0.0706 0.0667 0.0627 0.0588 0.0549 0.0510 0.0471 0.0431 0.0392 0.0353 0.0314 0.0275 0.0235 0.0196 0.0157 0.0118 0.0078 0.0039 0.0000];
        blue=[ 1.0000 0.9961 0.9922 0.9882 0.9843 0.9804 0.9765 0.9725 0.9686 0.9647 0.9608 0.9569 0.9529 0.9490 0.9451 0.9412 0.9373 0.9333 0.9294 0.9255 0.9216 0.9176 0.9137 0.9098 0.9059 0.9020 0.8980 0.8941 0.8902 0.8863 0.8824 0.8784 0.8745 0.8706 0.8667 0.8627 0.8588 0.8549 0.8510 0.8471 0.8431 0.8392 0.8353 0.8314 0.8275 0.8235 0.8196 0.8157 0.8118 0.8078 0.8039 0.8000 0.7961 0.7922 0.7882 0.7843 0.7804 0.7765 0.7725 0.7686 0.7647 0.7608 0.7569 0.7529 0.7490 0.7451 0.7412 0.7373 0.7333 0.7294 0.7255 0.7216 0.7176 0.7137 0.7098 0.7059 0.7020 0.6980 0.6941 0.6902 0.6863 0.6824 0.6784 0.6745 0.6706 0.6667 0.6627 0.6588 0.6549 0.6510 0.6471 0.6431 0.6392 0.6353 0.6314 0.6275 0.6235 0.6196 0.6157 0.6118 0.6078 0.6039 0.6000 0.5961 0.5922 0.5882 0.5843 0.5804 0.5765 0.5725 0.5686 0.5647 0.5608 0.5569 0.5529 0.5490 0.5451 0.5412 0.5373 0.5333 0.5294 0.5255 0.5216 0.5176 0.5137 0.5098 0.5059 0.5020 0.4980 0.4941 0.4902 0.4863 0.4824 0.4784 0.4745 0.4706 0.4667 0.4627 0.4588 0.4549 0.4510 0.4471 0.4431 0.4392 0.4353 0.4314 0.4275 0.4235 0.4196 0.4157 0.4118 0.4078 0.4039 0.4000 0.3961 0.3922 0.3882 0.3843 0.3804 0.3765 0.3725 0.3686 0.3647 0.3608 0.3569 0.3529 0.3490 0.3451 0.3412 0.3373 0.3333 0.3294 0.3255 0.3216 0.3176 0.3137 0.3098 0.3059 0.3020 0.2980 0.2941 0.2902 0.2863 0.2824 0.2784 0.2745 0.2706 0.2667 0.2627 0.2588 0.2549 0.2510 0.2471 0.2431 0.2392 0.2353 0.2314 0.2275 0.2235 0.2196 0.2157 0.2118 0.2078 0.2039 0.2000 0.1961 0.1922 0.1882 0.1843 0.1804 0.1765 0.1725 0.1686 0.1647 0.1608 0.1569 0.1529 0.1490 0.1451 0.1412 0.1373 0.1333 0.1294 0.1255 0.1216 0.1176 0.1137 0.1098 0.1059 0.1020 0.0980 0.0941 0.0902 0.0863 0.0824 0.0784 0.0745 0.0706 0.0667 0.0627 0.0588 0.0549 0.0510 0.0471 0.0431 0.0392 0.0353 0.0314 0.0275 0.0235 0.0196 0.0157 0.0118 0.0078 0.0039 0.0000];
    case 'greyscale'
        red=[ 0.0000 0.0039 0.0078 0.0118 0.0157 0.0196 0.0235 0.0275 0.0314 0.0353 0.0392 0.0431 0.0471 0.0510 0.0549 0.0588 0.0627 0.0667 0.0706 0.0745 0.0784 0.0824 0.0863 0.0902 0.0941 0.0980 0.1020 0.1059 0.1098 0.1137 0.1176 0.1216 0.1255 0.1294 0.1333 0.1373 0.1412 0.1451 0.1490 0.1529 0.1569 0.1608 0.1647 0.1686 0.1725 0.1765 0.1804 0.1843 0.1882 0.1922 0.1961 0.2000 0.2039 0.2078 0.2118 0.2157 0.2196 0.2235 0.2275 0.2314 0.2353 0.2392 0.2431 0.2471 0.2510 0.2549 0.2588 0.2627 0.2667 0.2706 0.2745 0.2784 0.2824 0.2863 0.2902 0.2941 0.2980 0.3020 0.3059 0.3098 0.3137 0.3176 0.3216 0.3255 0.3294 0.3333 0.3373 0.3412 0.3451 0.3490 0.3529 0.3569 0.3608 0.3647 0.3686 0.3725 0.3765 0.3804 0.3843 0.3882 0.3922 0.3961 0.4000 0.4039 0.4078 0.4118 0.4157 0.4196 0.4235 0.4275 0.4314 0.4353 0.4392 0.4431 0.4471 0.4510 0.4549 0.4588 0.4627 0.4667 0.4706 0.4745 0.4784 0.4824 0.4863 0.4902 0.4941 0.4980 0.5020 0.5059 0.5098 0.5137 0.5176 0.5216 0.5255 0.5294 0.5333 0.5373 0.5412 0.5451 0.5490 0.5529 0.5569 0.5608 0.5647 0.5686 0.5725 0.5765 0.5804 0.5843 0.5882 0.5922 0.5961 0.6000 0.6039 0.6078 0.6118 0.6157 0.6196 0.6235 0.6275 0.6314 0.6353 0.6392 0.6431 0.6471 0.6510 0.6549 0.6588 0.6627 0.6667 0.6706 0.6745 0.6784 0.6824 0.6863 0.6902 0.6941 0.6980 0.7020 0.7059 0.7098 0.7137 0.7176 0.7216 0.7255 0.7294 0.7333 0.7373 0.7412 0.7451 0.7490 0.7529 0.7569 0.7608 0.7647 0.7686 0.7725 0.7765 0.7804 0.7843 0.7882 0.7922 0.7961 0.8000 0.8039 0.8078 0.8118 0.8157 0.8196 0.8235 0.8275 0.8314 0.8353 0.8392 0.8431 0.8471 0.8510 0.8549 0.8588 0.8627 0.8667 0.8706 0.8745 0.8784 0.8824 0.8863 0.8902 0.8941 0.8980 0.9020 0.9059 0.9098 0.9137 0.9176 0.9216 0.9255 0.9294 0.9333 0.9373 0.9412 0.9451 0.9490 0.9529 0.9569 0.9608 0.9647 0.9686 0.9725 0.9765 0.9804 0.9843 0.9882 0.9922 0.9961 1.0000];
        green=[ 0.0000 0.0039 0.0078 0.0118 0.0157 0.0196 0.0235 0.0275 0.0314 0.0353 0.0392 0.0431 0.0471 0.0510 0.0549 0.0588 0.0627 0.0667 0.0706 0.0745 0.0784 0.0824 0.0863 0.0902 0.0941 0.0980 0.1020 0.1059 0.1098 0.1137 0.1176 0.1216 0.1255 0.1294 0.1333 0.1373 0.1412 0.1451 0.1490 0.1529 0.1569 0.1608 0.1647 0.1686 0.1725 0.1765 0.1804 0.1843 0.1882 0.1922 0.1961 0.2000 0.2039 0.2078 0.2118 0.2157 0.2196 0.2235 0.2275 0.2314 0.2353 0.2392 0.2431 0.2471 0.2510 0.2549 0.2588 0.2627 0.2667 0.2706 0.2745 0.2784 0.2824 0.2863 0.2902 0.2941 0.2980 0.3020 0.3059 0.3098 0.3137 0.3176 0.3216 0.3255 0.3294 0.3333 0.3373 0.3412 0.3451 0.3490 0.3529 0.3569 0.3608 0.3647 0.3686 0.3725 0.3765 0.3804 0.3843 0.3882 0.3922 0.3961 0.4000 0.4039 0.4078 0.4118 0.4157 0.4196 0.4235 0.4275 0.4314 0.4353 0.4392 0.4431 0.4471 0.4510 0.4549 0.4588 0.4627 0.4667 0.4706 0.4745 0.4784 0.4824 0.4863 0.4902 0.4941 0.4980 0.5020 0.5059 0.5098 0.5137 0.5176 0.5216 0.5255 0.5294 0.5333 0.5373 0.5412 0.5451 0.5490 0.5529 0.5569 0.5608 0.5647 0.5686 0.5725 0.5765 0.5804 0.5843 0.5882 0.5922 0.5961 0.6000 0.6039 0.6078 0.6118 0.6157 0.6196 0.6235 0.6275 0.6314 0.6353 0.6392 0.6431 0.6471 0.6510 0.6549 0.6588 0.6627 0.6667 0.6706 0.6745 0.6784 0.6824 0.6863 0.6902 0.6941 0.6980 0.7020 0.7059 0.7098 0.7137 0.7176 0.7216 0.7255 0.7294 0.7333 0.7373 0.7412 0.7451 0.7490 0.7529 0.7569 0.7608 0.7647 0.7686 0.7725 0.7765 0.7804 0.7843 0.7882 0.7922 0.7961 0.8000 0.8039 0.8078 0.8118 0.8157 0.8196 0.8235 0.8275 0.8314 0.8353 0.8392 0.8431 0.8471 0.8510 0.8549 0.8588 0.8627 0.8667 0.8706 0.8745 0.8784 0.8824 0.8863 0.8902 0.8941 0.8980 0.9020 0.9059 0.9098 0.9137 0.9176 0.9216 0.9255 0.9294 0.9333 0.9373 0.9412 0.9451 0.9490 0.9529 0.9569 0.9608 0.9647 0.9686 0.9725 0.9765 0.9804 0.9843 0.9882 0.9922 0.9961 1.0000];
        blue=[ 0.0000 0.0039 0.0078 0.0118 0.0157 0.0196 0.0235 0.0275 0.0314 0.0353 0.0392 0.0431 0.0471 0.0510 0.0549 0.0588 0.0627 0.0667 0.0706 0.0745 0.0784 0.0824 0.0863 0.0902 0.0941 0.0980 0.1020 0.1059 0.1098 0.1137 0.1176 0.1216 0.1255 0.1294 0.1333 0.1373 0.1412 0.1451 0.1490 0.1529 0.1569 0.1608 0.1647 0.1686 0.1725 0.1765 0.1804 0.1843 0.1882 0.1922 0.1961 0.2000 0.2039 0.2078 0.2118 0.2157 0.2196 0.2235 0.2275 0.2314 0.2353 0.2392 0.2431 0.2471 0.2510 0.2549 0.2588 0.2627 0.2667 0.2706 0.2745 0.2784 0.2824 0.2863 0.2902 0.2941 0.2980 0.3020 0.3059 0.3098 0.3137 0.3176 0.3216 0.3255 0.3294 0.3333 0.3373 0.3412 0.3451 0.3490 0.3529 0.3569 0.3608 0.3647 0.3686 0.3725 0.3765 0.3804 0.3843 0.3882 0.3922 0.3961 0.4000 0.4039 0.4078 0.4118 0.4157 0.4196 0.4235 0.4275 0.4314 0.4353 0.4392 0.4431 0.4471 0.4510 0.4549 0.4588 0.4627 0.4667 0.4706 0.4745 0.4784 0.4824 0.4863 0.4902 0.4941 0.4980 0.5020 0.5059 0.5098 0.5137 0.5176 0.5216 0.5255 0.5294 0.5333 0.5373 0.5412 0.5451 0.5490 0.5529 0.5569 0.5608 0.5647 0.5686 0.5725 0.5765 0.5804 0.5843 0.5882 0.5922 0.5961 0.6000 0.6039 0.6078 0.6118 0.6157 0.6196 0.6235 0.6275 0.6314 0.6353 0.6392 0.6431 0.6471 0.6510 0.6549 0.6588 0.6627 0.6667 0.6706 0.6745 0.6784 0.6824 0.6863 0.6902 0.6941 0.6980 0.7020 0.7059 0.7098 0.7137 0.7176 0.7216 0.7255 0.7294 0.7333 0.7373 0.7412 0.7451 0.7490 0.7529 0.7569 0.7608 0.7647 0.7686 0.7725 0.7765 0.7804 0.7843 0.7882 0.7922 0.7961 0.8000 0.8039 0.8078 0.8118 0.8157 0.8196 0.8235 0.8275 0.8314 0.8353 0.8392 0.8431 0.8471 0.8510 0.8549 0.8588 0.8627 0.8667 0.8706 0.8745 0.8784 0.8824 0.8863 0.8902 0.8941 0.8980 0.9020 0.9059 0.9098 0.9137 0.9176 0.9216 0.9255 0.9294 0.9333 0.9373 0.9412 0.9451 0.9490 0.9529 0.9569 0.9608 0.9647 0.9686 0.9725 0.9765 0.9804 0.9843 0.9882 0.9922 0.9961 1.0000];
        
    otherwise
        warning(['Colour scheme "' name '" not defined. Ignoring.']);
        red =[];
        green = [];
        blue = [];
        
end

output = [red' green' blue'];

end

function output = zlibdecode(input, outType)
%ZLIBDECODE Decompress input bytes using ZLIB.
%
%    output = zlibdecode(input, outType)
%
% The function takes a compressed byte array INPUT and returns inflated
% bytes OUTPUT. The INPUT is a zlib compressed (u)int8 array. The OUTPUT is
% always an 1-by-N outType array. If unspecified OUTTYPE is taken to be
% uint8. JAVA must be enabled to use the function.

narginchk(1, 2);

if nargin == 1
    outType = 'uint8';
end

error(javachk('jvm'));

if ischar(input)
    warning('zlibdecode:inputTypeMismatch', ...
        'Input is char, but treated as uint8.');
    input = uint8(input);
end

if ~isa(input, 'int8') && ~isa(input, 'uint8')
    error('Input must be either int8 or uint8.');
end

buffer = java.io.ByteArrayOutputStream(); % Stores an array of bytes

% Anything written to zlib is decompressed and written to buffer
zlib = java.util.zip.InflaterOutputStream(buffer);

zlib.write(input, 0, numel(input)); % Write the contents of input to zlib
                                    % I.e. decompress and store in buffer
zlib.close(); % Clean up

% Typecase from bytes to what we actually wanted
output = typecast(buffer.toByteArray(), outType)';

end

function output = makeValidFieldName(candidate,structure)
% Creates a unique and valid field name
narginchk(1,2);
if nargin < 2
    fNames = {};
else
    fNames = fieldnames(structure);
end

% Remove trailing spaces and nulls before processing
candidate = deblank(candidate);
% candidate = matlab.lang.makeValidName(candidate, ...
%     'Prefix', 'F_');
output = candidate;%matlab.lang.makeUniqueStrings(candidate, fNames);

end