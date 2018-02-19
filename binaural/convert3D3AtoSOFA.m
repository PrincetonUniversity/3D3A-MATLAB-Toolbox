function convert3D3AtoSOFA(matFile,SOFAFile,varargin)
%CONVERT3D3ATOSOFA Convert HRIRs stored in a MAT file in 3D3A format 
%   into a SOFA file.
%   CONVERT3D3ATOSOFA(matFile,SOFAFile) loads the HRIRs stored in matFile
%   in the 3D3A format and saves them to SOFAFile using the defaults
%   specified below.
%
%   CONVERT3D3ATOSOFA(...,Name1,Value1,Name2,Value2) specifies optional 
%   comma-separated pairs of Name,Value arguments, where Name is the 
%   argument name and Value is the corresponding value. Name must appear 
%   inside single quotes (' '). You can specify several name and value pair
%   arguments in any order as Name1,Value1,...,NameN,ValueN. Valid 
%   Name,Value arguments are as follows:
%   
%   'variables'     cell array of character vectors that specify the 
%                   variable names to be converted. Possible options 
%                   include (with option 1 being default):
%                   1. {'hrirL','hrirR','sourcePosition','sampleRate'}
%                   2. {'birL','birR','sourcePosition','sampleRate'}
%                   3. {'refL','refR','referencePosition','sampleRate'}
%
%   'coordinate'    character array that should either specify 'cartesian'
%                   (default) or 'spherical'. The matFile must have source 
%                   positions specified in cartesian coordinates, so this 
%                   option dictates the coordinate system to be used in the
%                   SOFAFile only.
%
%   'sample rate'   sampling rate in Hz. IRs in the matFile are resample to
%                   the specified sampling rate before exporting to the
%                   SOFAFile. The default is the sampling rate of the IRs
%                   in the matFile.
%
%   'mics'          string specifying the make and model of the binaural 
%                   microphones used. The default is 'Theoretica Applied 
%                   Physics BACCH-BM Pro 176-177'.
%
%   'sex'           string specifying the sex of subject. The default is
%                   'Unknown'.
%
%   'age'           string specifying the age of the subject. The default
%                   is 'Unknown'.
%
%   'date'          string specifying the measurement date. The default is
%                   the date of creation of the SOFAFile.

%   ==============================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2017 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

S = load(matFile);

[filePath, fileName, ~] = fileparts(matFile);
if nargin < 2 || isempty(SOFAFile)
    SOFAPath = filePath;
    SOFAName = fileName;
else
    [SOFAPath, SOFAName, ~] = fileparts(SOFAFile);
    if isempty(SOFAPath) || ~isdir(SOFAPath)
        SOFAPath = filePath;
        warning('Bad path for SOFA file; using MAT file path.');
    end
end

% Use alternate variables
indx = find(strcmpi(varargin,'variables'),1,'first');
if indx
    varNames = varargin{indx+1};
else
    varNames = {'hrirL','hrirR','sourcePosition','sampleRate'};
end
hrirL = S.(varNames{1});
hrirR = S.(varNames{2});
sourcePosition = S.(varNames{3});
sampleRate = S.(varNames{4});

% Specify coordinate system
indx = find(strcmpi(varargin,'coordinate'),1,'first');
if indx
    coord = varargin{indx+1};
else
    coord = 'cartesian';
end

% Resample if necessary
indx = find(strcmpi(varargin,'sample rate'),1,'first');
if indx
    newSampleRate = varargin{indx+1};
    if ~isequal(newSampleRate,sampleRate)
        hrirL = resample(hrirL,newSampleRate,sampleRate);
        hrirR = resample(hrirR,newSampleRate,sampleRate);
        sampleRate = newSampleRate;
    end
end

% Specify binaural microphones used
indx = find(strcmpi(varargin,'mics'),1,'first');
if indx
    micName = varargin{indx+1};
else
    micName = 'Theoretica Applied Physics BACCH-BM Pro 176-177';
end

% Specify subject sex
indx = find(strcmpi(varargin,'sex'),1,'first');
if indx
    sex = varargin{indx+1};
else
    sex = 'Unknown';
end

% Specify subject age
indx = find(strcmpi(varargin,'age'),1,'first');
if indx
    age = varargin{indx+1};
else
    age = 'Unknown';
end

% Get an empty conventions structure and populate convention-specific info.
if strcmpi(varNames{1},'hrirL')
    Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
    SOFAName = [SOFAName,'_HRIRs'];
    Obj.GLOBAL_Title = 'Measured and equalized HRIRs';
    % Listener attributes
    Obj.GLOBAL_ListenerShortName = SOFAName;
    Obj.GLOBAL_ListenerDescription = ['Sex: ',sex,'; Age: ',age];
elseif strcmpi(varNames{1},'birL')
    Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
    SOFAName = [SOFAName,'_BIRs'];
    Obj.GLOBAL_Title = 'Measured BIRs';
    % Listener attributes
    Obj.GLOBAL_ListenerShortName = SOFAName;
    Obj.GLOBAL_ListenerDescription = ['Sex: ',sex,'; Age: ',age];
else
    Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
    SOFAName = [SOFAName,'_RIRs'];
    Obj.GLOBAL_Title = 'Measured RIRs';
end
SOFAFile = fullfile(SOFAPath, [SOFAName,'.sofa']);

% Add the IRs
Obj.Data.IR = cat(2,shiftdim(shiftdim(hrirL,-1),2),...
    shiftdim(shiftdim(hrirR,-1),2));
Obj.Data.SamplingRate = sampleRate;

% Add the positons
emitterPosition = sourcePosition((round(sourcePosition(:,2),3) == 0) ...
    & (sourcePosition(:,1) >= 0),:);
switch lower(coord)
    case 'cartesian'
        Obj.SourcePosition = sourcePosition;
        Obj.SourcePosition_Type = 'cartesian';
        Obj.SourcePosition_Units = 'metre';
        Obj.EmitterPosition = emitterPosition;
        Obj.EmitterPosition_Type = 'cartesian';
        Obj.EmitterPosition_Units = 'metre';
    case 'spherical'
        Obj.SourcePosition = sofaC2sofaS(sourcePosition);
        Obj.SourcePosition_Type = 'spherical';
        Obj.SourcePosition_Units = 'degree, degree, metre';
        Obj.EmitterPosition = sofaC2sofaS(emitterPosition);
        Obj.EmitterPosition_Type = 'spherical';
        Obj.EmitterPosition_Units = 'degree, degree, metre';
end

% Global attributes
Obj.GLOBAL_ApplicationName = 'MATLAB';
Obj.GLOBAL_ApplicationVersion = version;
Obj.GLOBAL_AuthorContact = ['rahulram@princeton.edu; ',...
    'josephgt@princeton.edu; choueiri@princeton.edu'];
Obj.GLOBAL_Comment = '';
% Specify measurement date
indx = find(strcmpi(varargin,'date'),1,'first');
if ~isempty(indx) && ~isempty(varargin{indx+1})
    Obj.GLOBAL_DateCreated = varargin{indx+1};
end
% TODO - Think about updating the comment field.
Obj.GLOBAL_History = 'Measured in the 3D3A lab.';
Obj.GLOBAL_License = 'No license provided; ask the author for permission';
% TODO - Discuss adding license with Eddie.
Obj.GLOBAL_Organization = ['3D Audio and Applied Acoustics (3D3A) Lab,',...
    ' Princeton University'];
Obj.GLOBAL_Origin = ['Acoustic impulse responses measured with ',...
    'multiple exponential sine sweep method.'];
Obj.GLOBAL_DatabaseName = '3D3A Lab HRTF Database';

% Room attributes
Obj.GLOBAL_RoomShortName = '3D3A lab anechoic chamber';
Obj.GLOBAL_RoomDescription = ['The 3D3A lab is ',...
    'equipped with a highly reconfigurable anechoic chamber that was ',...
    'built in-house and that has adjustable/removable walls which can ',...
    'be moved on casters. The walls, floor, and ceiling are covered ',...
    'with 1-foot melamine square wedges. The wedges are 8 inches ',...
    'deep and give the chamber a cutoff frequency of 425 Hz. The ',...
    'chamber is lit with strings of green LEDs that give it a ',...
    'pleasant ambiance. The lighting and the fact that the chamber is ',...
    'not soundproof make it a pleasant environment for human test ',...
    'subjects.'];
Obj.GLOBAL_RoomLocation = 'Princeton, NJ, USA';
Obj = SOFAaddVariable(Obj,'RoomDimensions','IC',[3.6,2.35,2.55]);
Obj = SOFAaddVariable(Obj,'RoomDimensions_Type','IS',['length, width, ',...
    'height']);
Obj = SOFAaddVariable(Obj,'RoomDimensions_Units','IS','metre');
Obj = SOFAaddVariable(Obj,'RoomVolume','I',23.1);
Obj = SOFAaddVariable(Obj,'RoomVolume_Units','IS','cubic metre');
Obj.GLOBAL_RoomGeometry = 'shoebox';

Obj = SOFAaddVariable(Obj,'RoomCornerA','IC',[-1.4,...
    -Obj.RoomDimensions(2)/2,-1.46]);
Obj = SOFAaddVariable(Obj,'RoomCornerA_Type','IS','cartesian');
Obj = SOFAaddVariable(Obj,'RoomCornerA_Units','IS','metre');
Obj = SOFAaddVariable(Obj,'RoomCornerB','IC',Obj.RoomDimensions + ...
    Obj.RoomCornerA);
Obj = SOFAaddVariable(Obj,'RoomCornerB_Type','IS','cartesian');
Obj = SOFAaddVariable(Obj,'RoomCornerB_Units','IS','metre');

% Hardware descriptions
Obj.GLOBAL_SourceShortName = 'HRTF measurement arc';
Obj.GLOBAL_SourceDescription = ['Nine loudspeakers mounted on a ',...
    'vertical arc at fixed elevations from -57 to +75 degrees. ',...
    'Signals are sent to each loudspeaker from a  D/A converter ',...
    '(RME M-32 DA).'];
Obj.GLOBAL_ReceiverShortName = micName;
Obj.GLOBAL_ReceiverDescription = ['In-ear microphones (',micName,...
    ') connected to the microphone pre-amps ',...
    'of a digital audio interface (RME MADIface XT). Left and right ',...
    'ears correspond to the first and second receivers, respectively.'];
Obj.GLOBAL_EmitterShortName = 'Genelec 8010A';
Obj.GLOBAL_EmitterDescription = ['Small, active two-way ',...
    '(bi-amplified), non-coaxial loudspeakers.'];

% Save the SOFA file
Obj = SOFAsave(SOFAFile, Obj);

end

