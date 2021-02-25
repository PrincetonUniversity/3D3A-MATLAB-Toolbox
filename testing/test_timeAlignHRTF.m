%% Testing the removeHRTFLinearPhase function

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
%   Copyright (c) 2021 Princeton University
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

%% Specify RS-HRTF parameters

% Head parameters
head.a = 0.09;
head.eL = [90,0];
% Source parameters
rhoVec = [2;inf];
rVec = (head.a)*rhoVec;
rVecLen = length(rVec);
source.azVec = (90:30:270).';
source.elVec = 0;
% Computation parameters
dsp.fS = 48000;
dsp.T = 0.01;
dsp.method = {'sridharchoueiri2020',inf,4};
dsp.causalflag = false;
dsp.pow2flag = false;

%% Compute RS-HRTF with center normalization

dsp.normloc = 'center';

hrtfCell_center = cell(rVecLen,1);
for ii = 1:rVecLen
    source.r = rVec(ii);
    fprintf('Computing for r = %4.3f...\n',source.r);
    [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
    hrtfCell_center{ii,1} = fft(head.hrirL);
end

clearvars ii

%% Compute RS-HRTF with center_tdc normalization

dsp.normloc = 'center_tdc';

hrtfCell_center_tdc = cell(rVecLen,1);
for ii = 1:rVecLen
    source.r = rVec(ii);
    fprintf('Computing for r = %4.3f...\n',source.r);
    [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
    hrtfCell_center_tdc{ii,1} = fft(head.hrirL);
end

clearvars ii

%% Compute RS-HRTF with center_tdc_corr normalization

dsp.normloc = 'center_tdc_corr';

hrtfCell_center_tdc_corr = cell(rVecLen,1);
for ii = 1:rVecLen
    source.r = rVec(ii);
    fprintf('Computing for r = %4.3f...\n',source.r);
    [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
    hrtfCell_center_tdc_corr{ii,1} = fft(head.hrirL);
end

clearvars ii

%% Apply phase correction to center norm RS-HRTF to get center_tdc norm

fVec = getFreqVec(dsp.fS,dsp.irLen);
nyqIndx = ceil((dsp.irLen+1)/2);
hrtfCell_center_tdc_check = cell(rVecLen,1);
hrtfCell_center_tdc_corr_check = cell(rVecLen,1);
for ii = 1:rVecLen
    hrtfCell_center_tdc_check{ii,1} = timeAlignHRTF(...
        hrtfCell_center{ii,1}(1:nyqIndx,:),source.Positions,...
        fVec(1:nyqIndx),head.a,head.eL,rVec(ii),'SridharChoueiri2021a');
    hrtfCell_center_tdc_corr_check{ii,1} = timeAlignHRTF(...
        hrtfCell_center{ii,1}(1:nyqIndx,:),source.Positions,...
        fVec(1:nyqIndx),head.a,head.eL,rVec(ii),'SridharChoueiri2021b');
end

clearvars ii

%% Compute HRTF error

center_tdc_err = cell(rVecLen,1);
center_tdc_corr_err = cell(rVecLen,1);
for ii = 1:rVecLen
    center_tdc_err{ii,1} = hrtfCell_center_tdc{ii,1}(1:nyqIndx,:)-...
        hrtfCell_center_tdc_check{ii,1};
    center_tdc_corr_err{ii,1} = hrtfCell_center_tdc_corr{ii,1}...
        (1:nyqIndx,:)-hrtfCell_center_tdc_corr_check{ii,1};
end

% All the following displayed values should be very close to 0.
% Ignore Nyquist due to way in which RS-HRTF is calculated there (see
% computeSphereHRTF function for more information)
disp(max(max(real(center_tdc_err{1,1}(1:end-1,:)))))
disp(max(max(imag(center_tdc_err{1,1}(1:end-1,:)))))
disp(max(max(real(center_tdc_err{2,1}(1:end-1,:)))))
disp(max(max(imag(center_tdc_err{2,1}(1:end-1,:)))))

disp(max(max(real(center_tdc_corr_err{1,1}(1:end-1,:)))))
disp(max(max(imag(center_tdc_corr_err{1,1}(1:end-1,:)))))
disp(max(max(real(center_tdc_corr_err{2,1}(1:end-1,:)))))
disp(max(max(imag(center_tdc_corr_err{2,1}(1:end-1,:)))))

clearvars ii
