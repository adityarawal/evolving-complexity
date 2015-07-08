function [ I1, I2, points_knn dist_knn nx ny] = fastKraskovMI( X, Y, k, varargin )
%KraskovMI computes the Kraskov estimator for the mutual information.
%   1. Input: X, Y (n x 1) vectors
%             k: nearest neighbour
%             zeroFix (optional): fix the negative estimation to 0 (default
%                                 false);
%
%   2. Output: I1, I2: the 2 estimator of MI, I(1), I(2) (see Ref.)
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

%Divide the space into a 2-D square grid of size (k/N)^(1/2)
%Store values/pointers of X/Y values of the random variables
twoD_blocklen = (k/nObs)^(1/2);
twoD_numblocks = ceil((max(X)-min(X))/twoD_blocklen); 
twoD_blocklen = (max(X)-min(X))/twoD_numblocks; %Update based on the rounded off value of num blocks
twoD_grid = cell(twoD_numblocks, twoD_numblocks);
twoD_gridloc = zeros(nObs, 2) -1; 
twoDgrid_x_min = 1;
twoDgrid_x_max = twoD_numblocks;
twoDgrid_y_min = 1;
twoDgrid_y_max = twoD_numblocks;
points_knn = zeros(nObs, 2);
dist_knn = zeros(nObs, 1);

%Create 1-D grid for x and y variables and compute number of points
%in axis which are strictly less than dist_knn(i) apart from x(i)
% and y(i) respectively
oneD_slicelen = (1/nObs);
oneD_numslice = ceil((max(X)-min(X))/oneD_slicelen);
oneD_slicelen = (max(X)-min(X))/oneD_numslice; %Update based on the rounded off value of num blocks
oneD_x_grid = cell(oneD_numslice);
oneD_x_gridloc = zeros(nObs, 2) -1; 
oneD_y_grid = cell(oneD_numslice);
oneD_y_gridloc = zeros(nObs, 2) -1; 
oneDgrid_x_min = 1;
oneDgrid_x_max = oneD_numslice;
oneDgrid_y_min = 1;
oneDgrid_y_max = oneD_numslice;
nx = zeros(nObs, 1); %Number of x points that are strictly nearer to the point of interest than the knn  
ny = zeros(nObs, 1); %Number of y points that are strictly nearer to the point of interest than the knn

for i = 1:nObs
    if (X(i) == 0)
        X(i) = X(i) + 0.0000001; %Add small number to make sure block 1 is seleced
    end
    if (Y(i) == 0)
        Y(i) = Y(i) + 0.0000001; %Add small number to make sure block 1 is selected
    end
    %Create 2D grid
    block_x_loc = ceil(X(i)/twoD_blocklen); %Change ceil to floor for python/c++
    block_y_loc = ceil(Y(i)/twoD_blocklen); %Change ceil to floor for python/c++
    twoD_grid{block_x_loc, block_y_loc} = [twoD_grid{block_x_loc, block_y_loc}; [X(i) Y(i)]];
    twoD_gridloc(i,:) = [block_x_loc block_y_loc];  
    
    %Create 1D grid
    slice_x_loc = ceil(X(i)/oneD_slicelen); %Change ceil to floor for python/c++
    slice_y_loc = ceil(Y(i)/oneD_slicelen); %Change ceil to floor for python/c++
    oneD_x_grid{slice_x_loc} = [oneD_x_grid{slice_x_loc} X(i)];
    oneD_y_grid{slice_y_loc} = [oneD_y_grid{slice_y_loc} Y(i)];
    oneD_x_gridloc(i) = slice_x_loc;  
    oneD_y_gridloc(i) = slice_y_loc;  
end

%For each point, find the k-nearest neighbor by checking nearby grids
for i=1:nObs
    found_knn = false;
    found_nn = false; 
    ring_num = 0;
    block_x_loc = twoD_gridloc(i,1);
    block_y_loc = twoD_gridloc(i,2); 
    list_nn = [];
    list_nndist = [];

    %Find the distance to the k nearest neighbor
    while (found_knn == false)

        %Step 1 - Evaluate current ring blocks and find the nearest point
        [list_of_blocks] = calc_ring_blocks(block_x_loc, block_y_loc, ring_num, ...
            twoDgrid_x_max, twoDgrid_x_min, twoDgrid_y_max, twoDgrid_y_min);
    
        [temp_list_nn, temp_list_nndist, found_nn] = twoD_nearest_neighbor(list_of_blocks, X(i), Y(i), twoD_grid,k);
        
        if(found_nn)
            list_nn = [list_nn; temp_list_nn];
            list_nndist = [list_nndist; temp_list_nndist];
        end
        
        %If nearest point found, then validate by going to Step 2.
        if size(list_nn,1)>=k
            
            %Step 2 - Validate by going to the next ring.
            ring_num = ring_num + 1;
            [list_of_blocks] = calc_ring_blocks(block_x_loc, block_y_loc, ring_num, ...
                twoDgrid_x_max, twoDgrid_x_min, twoDgrid_y_max, twoDgrid_y_min);
            
            [temp_list_nn, temp_list_nndist, found_nn] = twoD_nearest_neighbor(list_of_blocks, ...
                X(i), Y(i), twoD_grid, k);
                      
            if(found_nn)
                list_nn = [list_nn; temp_list_nn];
                list_nndist = [list_nndist; temp_list_nndist];
            end
            [sorted_list_nndist, sorted_ind] = sort(list_nndist, 'ascend');
            points_knn(i,:) = list_nn(sorted_ind(k),:);
            dist_knn(i) = sorted_list_nndist(k);
            found_knn = true;
            
        %If no point found, then go to the next ring and go back to Step 1    
        else       
            ring_num = ring_num + 1;   
        end                
    end
    
    %Find points on x grid which are strictly less than k nearest distance 
    %apart from the point of interest
    nx(i) = total_nearby_points(X(i), oneD_x_grid, dist_knn(i), oneDgrid_x_min, oneDgrid_x_max, oneD_slicelen);
    
    %Find points on y grid ....
    ny(i) = total_nearby_points(Y(i), oneD_y_grid, dist_knn(i), oneDgrid_y_min, oneDgrid_y_max, oneD_slicelen);

end

% mutual information estimators
I1 = psi(k) - sum(psi(nx + 1) + psi(ny + 1)) / nObs + psi(nObs);
%I2 = psi(k) - 1/k - sum(psi(nx2) + psi(ny2)) / nObs + psi(nObs);
I2 = 0;
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
