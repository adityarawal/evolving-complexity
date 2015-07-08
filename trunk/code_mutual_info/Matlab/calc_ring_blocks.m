function [list_of_blocks] = calc_ring_blocks(block_x_loc, block_y_loc, ring_num, ...
    grid_x_max, grid_x_min, grid_y_max, grid_y_min)

% Returns x, y coordinate of all the blocks in a specific ring around the
% current block

list_of_blocks = [];

for i=-ring_num:ring_num
    for j=-ring_num:ring_num
        next_block_x = block_x_loc + i;
        next_block_y = block_y_loc + j;
        %%if (((i== -(ring_num) || i== ring_num || i == 0) && (j== -(ring_num) || j== ring_num || j==0)) ...
        %%        && (i~=0 || j ~= 0) ...
        if ((i~= -ring_num && i ~= ring_num) && (j~= -ring_num && j ~= ring_num) ...
                || next_block_x < grid_x_min || next_block_x > grid_x_max ...
                || next_block_y < grid_y_min || next_block_y > grid_y_max)
            continue;
        else
            list_of_blocks = [list_of_blocks; next_block_x next_block_y];
        end
    end
end


return;