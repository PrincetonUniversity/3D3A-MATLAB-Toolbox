function [] = start3D3AMATLABToolbox()
%START3D3AMATLABTOOLBOX Start the 3D3A MATLAB Toolbox.

%   =======================================================================
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
%   Copyright (c) 2018 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a
%   copy of this software and associated documentation files (the 
%   "Software"), to deal in the Software without restriction, including 
%   without limitation the rights to use, copy, modify, merge, publish, 
%   distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to 
%   the following conditions:
%   
%   The above copyright notice and this permission notice shall be included
%   in all copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%   =======================================================================

% First, try to find LTFAT
foundLTFAT = false;
if exist('ltfatstart','file') == 2
    foundLTFAT = true;
else
    ltfatDir = dir(fullfile(userpath,'**','ltfatstart.m'));
    if ~isempty(ltfatDir)
        addpath(ltfatDir(1).folder)
        foundLTFAT = true;
    else
        warning(['Could not find LTFAT! Please add LTFAT to the MATLAB'...
            'search path, otherwise some functions may not work.']);
    end
end

% Next, try to find SOFA
foundSOFA = false;
if exist('SOFAstart','file') == 2
    foundSOFA = true;
else
    sofaDir = dir(fullfile(userpath,'**','SOFAstart.m'));
    if ~isempty(sofaDir)
        addpath(sofaDir(1).folder)
        foundSOFA = true;
    else
        warning(['Could not find SOFA! Please add SOFA to the MATLAB'...
            'search path, otherwise some functions may not work.']);
    end
end

% Next, try to find SFS
foundSFS = false;
if exist('SFS_start','file') == 2
    foundSFS = true;
else
    sfsDir = dir(fullfile(userpath,'**','SFS_start.m'));
    if ~isempty(sfsDir)
        addpath(sfsDir(1).folder)
        foundSFS = true;
    else
        warning(['Could not find SFS! Please add SFS to the MATLAB'...
            'search path, otherwise some functions may not work.']);
    end
end

% Next, try to find AMT
foundAMT = false;
if exist('amt_start','file') == 2
    foundAMT = true;
else
    amtDir = dir(fullfile(userpath,'**','amt_start.m'));
    if ~isempty(amtDir)
        addpath(amtDir(1).folder)
        foundAMT = true;
    else
        warning(['Could not find AMT! Please add AMT to the MATLAB'...
            'search path, otherwise some functions may not work.']);
    end
end

if foundAMT
    amt_start; % This should also start SOFA, LTFAT, and SFS.
    if foundLTFAT
        % Delete rms.m from the LTFAT toolbox because of naming conflict
        ltfatPath = fileparts(which('ltfatstart.m'));
        if exist(fullfile(ltfatPath,'sigproc','rms.m'),'file')
            delete(fullfile(ltfatPath,'sigproc','rms.m'));
        end
    end
    if foundSFS
        % Delete rms.m from the SFS toolbox because of naming conflict
        sfsPath = fileparts(which('SFS_start.m'));
        if exist(fullfile(sfsPath,'SFS_general','rms.m'),'file')
            delete(fullfile(sfsPath,'SFS_general','rms.m'));
        end
    end
else
    if foundSOFA
        SOFAstart;
    end
    if foundLTFAT
        ltfatstart;
        % Delete rms.m from the LTFAT toolbox because of naming conflict
        ltfatPath = fileparts(which('ltfatstart.m'));
        if exist(fullfile(ltfatPath,'sigproc','rms.m'),'file')
            delete(fullfile(ltfatPath,'sigproc','rms.m'));
        end
    end
    if foundSFS
        SFS_start;
        % Delete rms.m from the SFS toolbox because of naming conflict
        sfsPath = fileparts(which('SFS_start.m'));
        if exist(fullfile(sfsPath,'SFS_general','rms.m'),'file')
            delete(fullfile(sfsPath,'SFS_general','rms.m'));
        end
    end
end

[toolboxDir,~,~] = fileparts(which('start3D3AMATLABToolbox'));
addpath(genpath(toolboxDir))
rmpath(genpath(fullfile(toolboxDir,'.git')))

disp('3D3A MATLAB Toolbox found and initialized.')

end
