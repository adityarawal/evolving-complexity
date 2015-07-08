function [ I1, I2, nn_points dist_nn nx1 ny1] = KraskovMI( X, Y, k, varargin )
%KraskovMI computes the Kraskov estimator for the mutual information.
%   1. Input: X, Y (n x 1) vectors
%             k: nearest neighbour
%             zeroFix (optional): fix the negative estimation to 0 (default
%                                 false);
%
%   2. Output: I1, I2: the two estimator of MI, I(1), I(2) (see Ref.)
%
% Ref: Kraskov, Alexander, Harald Stgbauer, and Peter Grassberger.
%      "Estimating mutual information." Physical review E 69.6 (2004): 066138.
%
% Author: Paolo Inglese <paolo.ingls@gmail.com>
% Last revision: 17-05-2015

if nargin < 3 || nargin > 4
    error('Wrong input number.');
end
if nargin == 3
    zeroFix = false;
end
if nargin == 4
    if ~islogical(varargin{1})
        error('zeroFix must be true or false');
    else
        zeroFix = varargin{1};
    end
end
    

if size(X, 1) ~= size(Y, 1)
    error('X and Y must contain the same number of samples');
end

nObs = size(X, 1);

% compute distance between each sample and its k-th nearest neighbour
dz = zeros(nObs, nObs);
dx = zeros(nObs, nObs);
dy = zeros(nObs, nObs);
for i = 1:nObs
    for j = 1:nObs
        dx(i,j) = pdist([X(i, :); X(j, :)]);
        dy(i,j) = pdist([Y(i, :); Y(j, :)]);
        dz(i,j) = max([dx(i, j), dy(i, j)]);
    end
end

% find nx(i) and ny(i)
Eps = zeros(nObs, 1);
Nn = zeros(nObs, 1);

nx1 = zeros(nObs, 1);
ny1 = zeros(nObs, 1);
nx2 = zeros(nObs, 1);
ny2 = zeros(nObs, 1);
nn_points = zeros(nObs, 2);
dist_points = zeros(nObs, 1);

for i = 1:nObs
    
    dxSample = dx(i, :);
    dxSample(i) = [];
    
    dySample = dy(i, :);
    dySample(i) = [];
    
    dzSample = dz(i, :);
    dzSample(i) = [];
    temp_X = X;
    temp_Y = Y;
    temp_X(i) = [];
    temp_Y(i) = [];
    [EpsSample, NnSample] = sort(dzSample, 'ascend');
    Eps(i) = EpsSample(k);
    Nn(i) = NnSample(k);
    nn_points(i, :) = [temp_X(Nn(i)) temp_Y(Nn(i))]; 
    dist_nn(i,1) = Eps(i);
    
    nx1(i) = sum(dxSample < Eps(i));
    ny1(i) = sum(dySample < Eps(i));
    
    nx2(i) = sum(dxSample <= Eps(i));
    ny2(i) = sum(dySample <= Eps(i));
    
end

% mutual information estimators
I1 = psi(k) - sum(psi(nx1 + 1) + psi(ny1 + 1)) / nObs + psi(nObs);
I2 = psi(k) - 1/k - sum(psi(nx2) + psi(ny2)) / nObs + psi(nObs);

if (zeroFix)
    if I1 < 0
        warning('First estimator is negative -> 0');
        I1 = 0;
    end
    if I2 < 0
        warning('Second estimator is negative -> 0');
        I2 = 0;
    end
end

end
