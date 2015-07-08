function [nx] = total_nearby_points(x, oneD_grid, dist_knn, oneDgrid_min, oneDgrid_max, oneD_slicelen)
%Find points on x which are strictly less than k nearest distance
%apart from the point of interest

ring_num = 0;
init_slice_loc = ceil(x/oneD_slicelen);
list_of_slices = [init_slice_loc];
found_left_side = true;
found_right_side = true;
flag = false;
nx = 0;

%No points can be nearest than the duplicate point
if (dist_knn == 0)
    nx = 0;
    return;
end 
    
while(found_left_side || found_right_side) %Check to ensure both sides of the point are explored
    
    for i=1:length(list_of_slices)
        slice_loc = list_of_slices(i);
        points_in_slice = oneD_grid{slice_loc};
        for j=1:length(points_in_slice)
            dx = pdist([points_in_slice(j); x]);
            if (dx == 0 && flag == false) %Ignore the point of interest but consider other duplicate points
                flag = true;
                continue;
            end
            if (dx < dist_knn)
                nx = nx + 1;
            else
                if (slice_loc < init_slice_loc) %Nothing found on left side
                    found_left_side = false;
                elseif (slice_loc > init_slice_loc) %Nothing found on right side
                    found_right_side = false;
                else
                    if (points_in_slice(j) < x)
                        found_left_side=false;
                    elseif (points_in_slice(j) > x)
                        found_right_side = false;
                    else
                        display('Error: Duplicate point in total_nearby_points()')
                    end                          
                end
             end
        end
    end
    
    ring_num = ring_num + 1;
    list_of_slices = [];
    if (((init_slice_loc-ring_num) >= oneDgrid_min) && found_left_side == true)
        list_of_slices = [list_of_slices (init_slice_loc-ring_num)];
    else
        found_left_side = false;
    end
    if (((init_slice_loc+ring_num) <= oneDgrid_max) && found_right_side == true)
        list_of_slices = [list_of_slices (init_slice_loc+ring_num)];
    else
        found_right_side = false;
    end
end

return;
