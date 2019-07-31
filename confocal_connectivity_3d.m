%Last modification: July 29th, 2019
%This function calculates the number of nodes, max/avg coordination number,
%length & radius, number of tails not on the image border and tail
%percentage from the input auto skeleton and spatial graph statistics files.
%Written by Geena. Modified by Edna

function [out] = confocal_connectivity_3d( skeleton_file, statistics_file)

num_nodes = xlsread(statistics_file,'Node Statistics','A2:A2'); %read total number of nodes
num_segments = xlsread(statistics_file,'Graph Summary','B2:B2');  %read total number of segments

%save data on node x,y coordinates and coordination number
nodes = xlsread(skeleton_file,1,['B2:E',num2str(num_nodes+1)]);

%save data on thickness at each point
point_thickness = xlsread(skeleton_file,2);
point_thickness = point_thickness(:,2);

%save data on segment length and mean radius
segments = xlsread(statistics_file,3,['B2:E',num2str(num_segments+1)]);
segments(:,2:3) = [];

%calculate max and average coordination numbers
coord = nodes(:,4);
max_coord = max(coord);
avg_coord = mean(coord(coord > 1));
stdv_coord=std(coord(coord > 1));

%calculate max and average segment lengths
max_length = max(segments(:,1));
avg_length = mean(segments(:,1));
stdv_length = std(segments(:,1));

%calculate max and average segment diameter
max_diam = max(point_thickness .* 2);
avg_diam = mean(point_thickness .* 2);
stdv_diam = std(point_thickness .* 2);

%count the number of tails not on the image border
tails_not_bb = 0;
for i = 1:num_nodes
    if nodes(i,4)== 1
        %count tail only if not on image border
        if nodes(i,1) == min(nodes(:,1)) || nodes(i,1) == max(nodes(:,1)) || ...
                nodes(i,2) == min(nodes(:,2)) || nodes(i,2) == max(nodes(:,2)) || ...
                nodes(i,3) == min(nodes(:,3)) || nodes(i,3) == max(nodes(:,3))
            continue
        end
        tails_not_bb = tails_not_bb + 1;
    end
end

%save data into return column vector
out = [num_nodes; max_coord; avg_coord; stdv_coord;max_length; avg_length; stdv_length; max_diam; ...
    avg_diam; stdv_diam; tails_not_bb; tails_not_bb/num_nodes];

end
