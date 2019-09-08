% test of “Five types of breast cancer” database

% What is sum of the following masks
what = 15;

% 1 - Initialise map
% 2 - Draw original distribution
% 4 - Map fitting for standard parameters
% 8 - Draw fitted map
% 16 - Draw fitted map in internal coordinates, nodes
% 32 - Draw fitted map in internal coordinates, edges
% 64 - Draw fitted map in internal coordinates, faces

% Classification in col vector
% Data in d1n matrix
% Colours to draw are defined now
colours = ['b'; 'g'];

% Data transformation if necessary
if exist('d1n', 'var') == 0
    %There is no data. Load
    load('breastCancer.mat');
end

if size(d1n, 1) > 286
    %we need to transpose matrix
    d1n = d1n';
end

if bitand(what, 1) ~= 0
    % Create and initialise map
    map = rect2DMap(30, 30);
    tic;
    init(map, d1n, 'pci');
    toc
end

if bitand(what, 2) ~= 0
    % Draw the firs mat to see distribution
    tic;
    drawMap(map, d1n, 'classes', col, 'markColour', colours,...
        'nodeMarker', 'none', 'lineWidth', 0.5);
    toc
end

if bitand(what, 4) ~= 0
    % Map fitting for standard parameters
    tic;
    EM(map, d1n, 'stretch', 0, 'bend', 0.1);
    toc
end

if bitand(what, 8) ~= 0
    % Draw the firs mat to see distribution
    tic;
    drawMap(map, d1n, 'classes', col, 'markColour', colours,...
        'nodeMarker', 'none', 'lineWidth', 0.5);
    toc
end

if bitand(what, 16) ~= 0
    % Draw fitted map in internal coordinates, nodes
%     for k=1:2
%         if k == 1
%             map = map1;
%         else
%             map = map2;
%         end
        tic;
        drawMapInt(map, d1n, 0, 'lineWidth', 0.5, 'nodeMarker', 'h',...
             'classes', col, 'markColour', colours);
        toc
%     end
end

if bitand(what, 32) ~= 0
    % Draw fitted map in internal coordinates, edges
    for k=1:2
        if k == 1
            map = map1;
        else
            map = map2;
        end
        tic;
        drawMapInt(map, d1n, 1, 'nodeMarker', 'none', 'lineWidth', 0.5,...
            'classes', col, 'markColour', colours,...
                'nodeMarker', 'none', 'lineWidth', 0.5);
        toc
    end
end

if bitand(what, 64) ~= 0
    % Draw fitted map in internal coordinates, faces
    for k=1:2
        if k == 1
            map = map1;
        else
            map = map2;
        end
        tic;
        drawMapInt(map, d1n, 2, 'nodeMarker', 'none', 'lineWidth', 0.5,...
            'classes', col, 'markColour', colours,...
                'nodeMarker', 'none', 'lineWidth', 0.5);
        toc
    end
end
