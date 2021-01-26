function combinedCell = csv2cell(filename)
%CSV2CELL Read a comma separated value file.
%   C = CSV2CELL(FILENAME) reads a comma separated value formatted file
%   FILENAME. The result is returned in the cell array C. C will be an M-by
%   -N cell array, where M is the number of lines in the CSV file and N is
%   the greatest number of comma-separated values in any single line.
%
%   See also CSVREAD.

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

fid = fopen(filename);
TLINE = fgetl(fid);
topLine = textscan(TLINE,'%s','delimiter',',');
headerRow = topLine{1}.';
nCols = length(headerRow);
csvFormat = strrep(strrep(int2str(ones(1,nCols)),' ',''),'1','%s');
notesData = textscan(fid,csvFormat,'delimiter',',');
fclose(fid);

nRows = length(notesData{1}) + 1;
for ii = 2:nCols
    temp = length(notesData{ii}) + 1;
    if temp < nRows
        for jj = 1:(nRows-temp)
            notesData{ii}{temp+jj-1} = '';
        end
%         nRows = temp;
    end
end
combinedCell = cell(nRows,nCols);
combinedCell(1,:) = headerRow;
for ii = 1:nCols
    combinedCell(2:nRows,ii) = notesData{ii}(1:nRows-1);
end

end
