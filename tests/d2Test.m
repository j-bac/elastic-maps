% test of �Five types of breast cancer� database

% What is sum of the masks below
what = 1;

% 1 - Initialise map
% 2 - Draw original distribution
% 4 - Map fitting for standard parameters
% 8 - Draw fitted map
% 16 - Draw fitted map in internal coordinates, nodes
% 32 - Draw fitted map in internal coordinates, edges
% 64 - Draw fitted map in internal coordinates, faces

% Classification in col vector
% Data in d2fn matrix
% Colours to draw are defined now
colours = ['b'; 'g'; 'm'];

% Data transformation if necessary
if exist('d2fn', 'var') == 0
    %There is no data. Load
    load('bladderCancer.mat');
end

if size(d2fn, 1) > 40
    %we need to transpose matrix
    d2fn = d2fn';
end

if bitand(what, 1) ~= 0
    % Create and initialise map
    map = rect2DMap(10, 10);
    tic;
    init(map, d2fn, 'pci');
    toc
end

if bitand(what, 1) ~= 0
    % Draw the firs mat to see distribution
    tic;
    drawMap(map, d2fn, 'classes', col, 'markColour', colours,...
        'nodeMarker', 'none', 'lineWidth', 0.5);
    toc
end

if bitand(what, 1) ~= 0
    % Map fitting for standard parameters
    tic;
    %EM(map, d2fn, 'stretch', 0, 'bend', 0.1);
    EM(map, d2fn, 'stretch', 0, 'bend', 0.1, 'potential', @L2,...
            'Number_of_intervals', 2, 'intshrinkage', 0.75);
    toc
end

if bitand(what, 1) ~= 0
    % Draw the fitted map to see distribution
    tic;
    drawMap(map, d2fn, 'classes', col, 'markColour', colours,...
        'nodeMarker', 'none', 'lineWidth', 0.5);
    toc
end
% 
% if bitand(what, 1) ~= 0
%     % Draw fitted map in internal coordinates, nodes
% %     for k=1:2
% %         if k == 1
% %             map = map1;
% %         else
% %             map = map2;
% %         end
%         tic;
%         drawMapInt(map, d2fn, 0, 'lineWidth', 0.5, 'nodeMarker', 'h',...
%              'classes', col, 'markColour', colours);
%         toc
% %     end
% end
% 
% if bitand(what, 1) ~= 0
%     % Draw fitted map in internal coordinates, edges
% %     for k=1:2
% %         if k == 1
% %             map = map1;
% %         else
% %             map = map2;
% %         end
%         tic;
%         drawMapInt(map, d2fn, 1, 'nodeMarker', 'none', 'lineWidth', 0.5,...
%             'classes', col, 'markColour', colours,...
%                 'nodeMarker', 'none', 'lineWidth', 0.5);
%         toc
% %     end
% end
% 
% if bitand(what, 1) ~= 0
%     % Draw fitted map in internal coordinates, faces
% %     for k=1:2
% %         if k == 1
% %             map = map1;
% %         else
% %             map = map2;
% %         end
%         tic;
%         drawMapInt(map, d2fn, 2, 'nodeMarker', 'none', 'lineWidth', 0.5,...
%             'classes', col, 'markColour', colours,...
%                 'nodeMarker', 'none', 'lineWidth', 0.5);
%         toc
% %     end
% end
