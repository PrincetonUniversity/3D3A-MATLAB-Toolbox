function pathOut = findFiles(pathIn,fileName)
%FINDFILES Recursively search for particular files in given path.
%   PO = FINDFILES(PI,F) recursively searches PI for files with name given
%   in F. If N such files are found, N > 0, the paths to the files are 
%   returned in the cell array PO. If no such files are found, PO is 
%   returned as an empty cell array.
%
%   See also FINDFOLDERS.

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

narginchk(2,2);

% Validate inputs
validateattributes(pathIn,{'char'},{'scalartext','nonempty'},...
    'findFiles','PI',1);
validateattributes(fileName,{'char'},{'scalartext','nonempty'},...
    'findFiles','F',2);

if exist(pathIn,'dir') ~= 7
    error('Input path, PI, not found.')
end

pathOut = {};
% Look for F in current directory
pathIn_foundItems = dir(fullfile(pathIn,fileName));
if ~isempty(pathIn_foundItems)
    % Ignore folders that may have the specified file name
    pathIn_foundFiles = pathIn_foundItems(~[pathIn_foundItems.isdir].');
    pathIn_numFoundFiles = size(pathIn_foundFiles,1);
    for ii = 1:pathIn_numFoundFiles
        pathToAppend = fullfile(pathIn_foundFiles(ii).folder,...
            pathIn_foundFiles(ii).name);
        pathOut = [pathOut; pathToAppend];
    end
end

% Next, look for F recursively in all sub-directories
% Begin by obtaining a listing of directories only
pathIn_contents = dir(pathIn); % Listing of all contents in input path
pathIn_dirsOnly = pathIn_contents([pathIn_contents.isdir].');
% Begin recursively searching directories for F
pathIn_size = size(pathIn_dirsOnly,1);
% If pathIn_size == 2, only '.' and '..' directories exist and can be
% ignored (it is assumed that '.' and '..' always exist)
if pathIn_size > 2
    % Loop through all folders in current directory
    for ii = 1:pathIn_size
        currentName = pathIn_dirsOnly(ii).name;
        % Ignore the '.' and '..' directories to prevent an infinite loop
        if ~strcmpi(currentName,'.') && ~strcmpi(currentName,'..')
            currentPath = fullfile(pathIn_dirsOnly(ii).folder,currentName);
            cellToAppend = findFiles(currentPath,fileName);
            pathOut = [pathOut; cellToAppend];
        end
    end
end

end
