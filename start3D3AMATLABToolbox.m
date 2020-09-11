function start3D3AMATLABToolbox()
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

% Check input count
narginchk(0,0)

disp('Starting 3D3A MATLAB Toolbox...')
[toolboxDir,~,~] = fileparts(which('start3D3AMATLABToolbox'));
depDir = fullfile(toolboxDir,'dependencies');

% First, try to find LTFAT
foundLTFAT = false;
if exist('ltfatstart','file') == 2
    foundLTFAT = true;
else
    ltfatDir1 = dir(fullfile(userpath,'**','ltfatstart.m'));
    ltfatDir2 = dir(fullfile(depDir,'**','ltfatstart.m'));
    if ~isempty(ltfatDir1)
        addpath(ltfatDir1(1).folder)
        foundLTFAT = true;
    elseif ~isempty(ltfatDir2)
        addpath(ltfatDir2(1).folder)
        foundLTFAT = true;
    else
        disp('Could not find local copy of dependency: LTFAT.')
        disp('Attempting to download from web...')
        if ispc
            if verLessThan('matlab','9.4')
                url = ['https://github.com/ltfat/ltfat/releases/',...
                    'download/2.4.0/ltfat-2.4.0-bin-win64.zip'];            
            else
                url = ['https://github.com/ltfat/ltfat/releases/',...
                    'download/2.4.0/ltfat-2.4.0-bin-win64-R2018a.zip'];
            end
            
            try
                unzip(url,fullfile(depDir,'ltfat'));
                foundLTFAT = true;
            catch
                disp('Unable to download/unzip LTFAT toolbox.')
            end
        else
            url = ['https://github.com/ltfat/ltfat/releases/download/',...
                '2.4.0/ltfat-2.4.0-src.tgz'];
            try
                gunzip(url,depDir);
                filename = 'ltfat-2.4.0-src';
                untar(fullfile(depDir,filename),fullfile(depDir,'ltfat'));
                delete(fullfile(depDir,filename))
                foundLTFAT = true;
            catch
                disp('Unable to download/unzip LTFAT toolbox.')
            end
        end
        
        if foundLTFAT
            ltfatDir = dir(fullfile(depDir,'**','ltfatstart.m'));
            if ~isempty(ltfatDir)
                addpath(ltfatDir(1).folder)
                disp(['LTFAT downloaded successfully and added to',...
                    ' search path.'])
            end
        end
    end
end

if ~foundLTFAT
    warning(['Could not find LTFAT! Please add LTFAT to the MATLAB',...
        ' search path, otherwise some functions may not work.']);
end

% Next, try to find SOFA
foundSOFA = false;
if exist('SOFAstart','file') == 2
    foundSOFA = true;
else
    sofaDir1 = dir(fullfile(userpath,'**','SOFAstart.m'));
    sofaDir2 = dir(fullfile(depDir,'**','SOFAstart.m'));
    if ~isempty(sofaDir1)
        addpath(sofaDir1(1).folder)
        foundSOFA = true;
    elseif ~isempty(sofaDir2)
        addpath(sofaDir2(1).folder)
        foundSOFA = true;    
    else
        disp('Could not find local copy of dependency: SOFA API.')
        disp('Attempting to download from web...')
        url = ['https://downloads.sourceforge.net/project/',...
            'sofacoustics/SOFA%20API%20for%20Matlab%20and%20Octave%',...
            '201.0.4.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%',...
            '2Fsofacoustics%2Ffiles%2FSOFA%2520API%2520for%2520Matlab',...
            '%2520and%2520Octave%25201.0.4.zip%2Fdownload&ts=1599856804'];
        zipfilepath = fullfile(depDir,'sofa.zip');
        try
            outfilename = websave(zipfilepath,url);
            try
                unzip(outfilename,fullfile(depDir,'sofa'));
                delete(outfilename)
                foundSOFA = true;
            catch
                disp('Unable to unzip downloaded SOFA API.')
            end
        catch
            disp('Unable to download SOFA API.')
        end
        
        if foundSOFA
            sofaDir = dir(fullfile(depDir,'**','SOFAstart.m'));
            if ~isempty(sofaDir)
                addpath(sofaDir(1).folder)
                disp(['SOFA API downloaded successfully and added to',...
                    ' search path.'])
            end
        end
    end
end

if ~foundSOFA
    warning(['Could not find the SOFA API! Please add the API to the',...
        ' MATLAB search path, otherwise some functions may not work.']);
end

% Next, try to find SFS
foundSFS = false;
if exist('SFS_start','file') == 2
    foundSFS = true;
else
    sfsDir1 = dir(fullfile(userpath,'**','SFS_start.m'));
    sfsDir2 = dir(fullfile(depDir,'**','SFS_start.m'));
    if ~isempty(sfsDir1)
        addpath(sfsDir1(1).folder)
        foundSFS = true;
    elseif ~isempty(sfsDir2)
        addpath(sfsDir2(1).folder)
        foundSFS = true; 
    else
        disp('Could not find local copy of the dependency: SFS Toolbox.')
        disp('Attempting to download from web...')
        url = ['https://github.com/sfstoolbox/sfs-matlab/archive/',...
            '2.5.0.zip'];
        try
            unzip(url,fullfile(depDir,'sfs'));
            foundSFS = true;
        catch
            disp('Unable to download/unzip SFS toolbox.')
        end
        
        if foundSFS
            sfsDir = dir(fullfile(depDir,'**','SFS_start.m'));
            if ~isempty(sfsDir)
                addpath(sfsDir(1).folder)
                disp(['SFS toolbox downloaded successfully and added',...
                    ' to search path.'])
            end
        end
    end
end

if ~foundSFS
    warning(['Could not find the SFS toolbox! Please add the toolbox',...
        ' to the MATLAB search path, otherwise some functions may not',...
        ' work.']);
end

% Next, try to find AMT
foundAMT = false;
if exist('amt_start','file') == 2
    foundAMT = true;
else
    amtDir1 = dir(fullfile(userpath,'**','amt_start.m'));
    amtDir2 = dir(fullfile(depDir,'**','amt_start.m'));
    if ~isempty(amtDir1)
        addpath(amtDir1(1).folder)
        foundAMT = true;
    elseif ~isempty(amtDir2)
        addpath(amtDir2(1).folder)
        foundAMT = true;
    else
        disp('Could not find local copy of the dependency: AMT.')
        disp('Attempting to download from web...')
        url = ['https://downloads.sourceforge.net/project/amtoolbox/',...
            'amtoolbox-full-0.10.0.zip?r=https%3A%2F%2Fsourceforge.net',...
            '%2Fprojects%2Famtoolbox%2Ffiles%2Flatest%2Fdownload%3F',...
            'source%3Dfiles&ts=1599858933'];
        zipfilepath = fullfile(depDir,'amt.zip');
        try
            outfilename = websave(zipfilepath,url);
            try
                unzip(outfilename,fullfile(depDir,'amt'));
                delete(outfilename)
                foundAMT = true;
            catch
                disp('Unable to unzip downloaded AMT.')
            end
        catch
            disp('Unable to download AMT.')
        end
        
        if foundAMT
            amtDir = dir(fullfile(depDir,'**','amt_start.m'));
            if ~isempty(amtDir)
                addpath(amtDir(1).folder)
                disp(['AMT downloaded successfully and added to search',...
                    ' path.'])
            end
        end
    end
end

if ~foundAMT
    warning(['Could not find the AMT! Please add the AMT to the MATLAB',...
        ' search path, otherwise some functions may not work.']);
else
    amtDir = which('amt_start.m');
    [amtDir,~,~] = fileparts(amtDir);
    amt_tpDir = dir(fullfile(amtDir,'**','thirdparty'));
    amt_tpDir = amt_tpDir(1).folder;
    if exist(fullfile(amt_tpDir,'ltfat'),'dir') == 7
        rmdir(fullfile(amt_tpDir,'ltfat'),'s');
    end
    
    if exist(fullfile(amt_tpDir,'sfs'),'dir') == 7
        rmdir(fullfile(amt_tpDir,'sfs'),'s');
    end
    
    if exist(fullfile(amt_tpDir,'SOFA'),'dir') == 7
        rmdir(fullfile(amt_tpDir,'SOFA'),'s');
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

addpath(genpath(toolboxDir))
rmpath(genpath(fullfile(toolboxDir,'.git')))

disp('3D3A MATLAB Toolbox initialization complete.')

end
