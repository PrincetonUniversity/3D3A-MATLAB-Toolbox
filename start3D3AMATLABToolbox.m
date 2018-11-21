function [] = start3D3AMATLABToolbox()
%START3D3AMATLABTOOLBOX Start the 3D3A MATLAB Toolbox.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
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

% Start the found toolbox(es)
if foundLTFAT && foundAMT
    amt_start; % Also starts LTFAT
elseif foundLTFAT && ~foundAMT
    ltfatstart;
elseif ~foundLTFAT && foundAMT
    try
        amt_start; % Might find LTFAT in AMT's thirdparty folder
    catch ME
        disp('Found AMT but could not start due to the following error:')
        messageLines = splitlines(ME.message);
        indentMessage = compose('\t %s',messageLines);
        disp(strjoin(indentMessage,'\n'))
    end
end

[toolboxDir,~,~] = fileparts(which('start3D3AMATLABToolbox'));
addpath(genpath(toolboxDir))
rmpath(genpath(fullfile(toolboxDir,'.git')))

disp('3D3A MATLAB Toolbox found and initialized.')

end