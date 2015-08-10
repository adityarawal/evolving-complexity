function [list_nn, list_nndist, found_nn] = twoD_nearest_neighbor(list_of_blocks, x, y, twoD_grid, k)

num_blocks = size(list_of_blocks,1);
% compute distance between each sample and its k-th nearest neighbour
% x, y - point of interest
% list_of_blocks - find the shortest distance between point of interest and
% the points residing inside the list_of_blocks

list_nndist = [];
list_nn = [];

%If the point of interest lies within the list_of_blocks, then exclude the point (not the
%duplicates) while calculating the nearest neighbor
flag = false;
for i = 1:size(list_of_blocks,1)
    block_x_loc = list_of_blocks(i,1);
    block_y_loc = list_of_blocks(i,2);
    points_in_block = twoD_grid{block_x_loc, block_y_loc};
    for j=1:size(points_in_block,1)
        dx = pdist([points_in_block(j, 1); x]);
        dy = pdist([points_in_block(j, 2); y]);
        if (dx == 0 && dy==0 && flag == false)
            flag = true;
            continue;
        else
            list_nndist = [list_nndist; max([dx, dy])];
            list_nn = [list_nn; [points_in_block(j, 1) points_in_block(j, 2)]];
        end
    end
end
    
if isempty(list_nn) %No point found
    found_nn = false;
    list_nn = [];
    list_nndist = [];
else %Some points found
    found_nn = true;
end

return