function [A,B,C,D,E] = ir2ss(h,N,METHOD)
%IR2SS Partial state-space realization of a finite impulse response.
%   [A,B,C,D,E] = IR2SS(H) returns a partial state-space realization of the 
%   finite impulse response (FIR) H using Beliczynski's realization 
%   algorithm [1].
%       For a single-input, single-output system, H may be a 1-by-L or
%       L-by-1 vector.
%       For a multi-input and/or multi-output system, H must be an array of
%       dimensions P-by-M-by-L, where P corresponds to the number of
%       outputs, M to the number of inputs, and L to the length of each
%       FIR.
%   By default, the order of the realization, N, is floor(L/2). D is always
%   returned as 0. Consequently, when computing an IR from the state-space
%   model, the first sample should be dropped. E is a structure containing
%   additional variables specific to the chosen METHOD (see below).
%   
%   [A,B,C,D,E] = IR2SS(H,N) optionally specifies the order to use. N must 
%   be a positive integer less than L. If N = L-1, the resulting state-
%   space realization will be in observable-canonical form. By default, the 
%   value of N is floor(L/2).
%
%   [A,B,C,D,E] = IR2SS(H,N,METHOD) optionally specifies the realization
%   algorithm to use. METHOD can take the following values: 'beliczynski'
%   (default), and 'silverman'.
%
%   See also SS2IR.

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

% Refs:
%   [1]. Antoulas (2005) - Ch. 4. Linear Dynamical Systems: Part 1 in 
%        Approximation of Large-Scale Dynamical Systems.

narginchk(1,3);

if isvector(h)
    hLen = length(h);
    bh = zeros(1,1,hLen);
    bh(1,1,:) = shiftdim(h);
    p = 1;
    m = 1;
elseif ndims(h) == 3
    [p,m,hLen] = size(h);
    bh = h;
else
    error('Invalid input dimensions.')
end

if nargin < 3 
    METHOD = 'beliczynski';
end

if nargin < 2 || isempty(N)
    N = floor(hLen/2);
end

validateattributes(N,{'double'},{'scalar','nonempty','nonnan','finite',...
    'integer','positive','<',hLen},'ir2ss','N',2);

% Generate partially-defined Hankel matrix
H = getBlockHankelMatrix(bh);
n = N*min([p,m]);
E = struct;
switch lower(METHOD)
    case 'silverman'
        rH = rank(H);
        if rH < n
            warning(['Rank of H is less than specified order. Reducing ',...
                'order to %d.'],rH);
            n = rH;
        end
        
        % Implement Silverman's realization algorithm
        phi = H(1:n,1:n);
        sPhi = H(1:n,(m+1):(m+n));
        Gamma = H(1:n,1:m);
        Lambda = H(1:p,1:n);
        
        A = phi\sPhi;
        B = phi\Gamma;
        C = Lambda;
    case 'beliczynski'
        [~,S,V] = svd(H);
        chLen = length(h(:));
        E.hsvs = diag(S);
        if n < chLen
            E.totbnd = 2*sum(E.hsvs(n+1:end));
        else
            E.totbnd = 0;
        end
        A = V(2:chLen,1:n)'*V(1:chLen-1,1:n);
        B = V(1,1:n)';
        C = (h(:)).'*V(1:chLen,1:n);
    otherwise
        error('Invalid METHOD specification.')
end

D = 0;

end
