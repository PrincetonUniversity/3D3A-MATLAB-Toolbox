function absPath = rel2abs(relPath)
%REL2ABS Convert relative path to absolute path.
%   B = REL2ABS(A) converts the relative path specified in A to an absolute 
%   path, B, where the absolute path is generated assuming that the
%   relative path is specified relative to the current folder in MATLAB. If
%   the path contains a shortcut for the user directory, this is also
%   replaced by the full path. The following shortcuts are recognized:
%       1. '~' on Unix-based systems (including macOS)
%       2. %userprofile% on Windows
%
%   See also GETUSERPATH.

%   =======================================================================
%   This file is part of the 3D3A MATLAB Toolbox.
%   
%   Contributing author(s), listed alphabetically by last name:
%   Rahulram Sridhar <rahulram@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2019 Princeton University
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

narginchk(1,1);

validateattributes(relPath,{'char'},{'scalartext','nonempty'},...
    'rel2abs','A',1);

cDir = pwd; % Current directory

cDir_sepIndxs = strfind(cDir,filesep); % File separator locations for cDir
relPath_DDotIndxs = strfind(relPath,'..');
relPath_numDDotIndxs = length(relPath_DDotIndxs);
if relPath_numDDotIndxs > 0
    absPath = [cDir(1:cDir_sepIndxs(end-relPath_numDDotIndxs+1)),...
        relPath(relPath_DDotIndxs(end)+2+length(filesep):end)];
else
    absPath = relPath;
end

if ismac || isunix
    userPathSym = '~';
elseif ispc
    userPathSym = '%userprofile%';
else
    disp('Platform not supported.')
end

userPathSymFlag = strfind(absPath,userPathSym);
if userPathSymFlag
    absPath = strrep(absPath,userPathSym,getUserPath);
end

end
