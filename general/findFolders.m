function pathOut = findFolders(pathIn,folderName,varargin)
%FINDFOLDERS Recursively search for particular folders in given path.
%   PO = FINDFOLDERS(PI,F) recursively searches PI for folders with name
%   given in F. If N such folders are found, N > 0, the paths to the 
%   folders are returned in the cell array PO. If no such folders are
%   found, PO is returned as an empty cell array. By default, FINDFOLDERS
%   returns all paths within the directory PI that contain a folder with
%   name given in F. For example, if searching for 'B' within the folder
%   structure: A/B/C/B, PO will contain the following directories:
%   1) A/B 2) A/B/C/B
%
%   PO = FINDFOLDERS(PI,F,'first') optionally stops the search as soon as
%   the first folder with name given in F is found. If another folder with
%   name F is located at any level within this folder, it will not be found
%   when this option is specified. For example, if searching for 'B' within 
%   the folder structure: A/B/C/B, PO will contain the following 
%   directories: 1) A/B
%
%   See also FINDFILES.

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
%   Copyright (c) 2020 Princeton University
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

narginchk(2,3);

% Validate inputs
validateattributes(pathIn,{'char'},{'scalartext','nonempty'},...
    'findFolders','PI',1);
validateattributes(folderName,{'char'},{'scalartext','nonempty'},...
    'findFolders','F',2);

indx = find(strcmpi(varargin,'first'),1);
if ~isempty(indx)
    firstFlag = true;
    firstString = 'first';
else
    firstFlag = false;
    firstString = 'default';
end

if exist(pathIn,'dir') ~= 7
    error('Input path, PI, not found.')
end

pathOut = {};
% Get listing of contents in input path
pathIn_contents = dir(pathIn);
% Extract listing of directories only
pathIn_dirsOnly = pathIn_contents([pathIn_contents.isdir].');
pathIn_size = size(pathIn_dirsOnly,1);
% If pathIn_size == 2, only '.' and '..' directories exist and can be
% ignored (it is assumed that '.' and '..' always exist)
if pathIn_size > 2
    for ii = 1:pathIn_size
        currentName = pathIn_dirsOnly(ii).name;
        currentDir = pathIn_dirsOnly(ii).folder;
        currentPath = fullfile(currentDir,currentName);
        % Ignore the '.' and '..' directories to prevent an infinite loop
        if ~strcmpi(currentName,'.') && ~strcmpi(currentName,'..')
            if strcmpi(currentName,folderName)
                pathOut = [pathOut;{currentPath}];
                if ~firstFlag
                    cellToAppend = findFolders(currentPath,folderName);
                    pathOut = [pathOut;cellToAppend];
                end
            else
                cellToAppend = findFolders(currentPath,folderName,...
                    firstString);
                pathOut = [pathOut;cellToAppend];
            end
        end
    end
end

end
