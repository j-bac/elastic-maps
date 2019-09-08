% test of “Five types of breast cancer” database in transpose shape

% Data transformation if necessary
if exist('d1n', 'var') == 0
    %There is no data. Load
    load('brestCancer.mat');
end

if size(d1n, 2) > 286
    %we need to transpose matrix
    d1n = d1n';
end

% Create and initialise map
fprintf('creation');
map = rect2DMap(30, 30);
tic;
init(map, d1n, 'pci');
toc

fprintf('Figure 1');
tic;
    drawMap(map, d1n,'nodeMarker', 'none', 'lineWidth', 0.5);
toc

    % Map fitting for standard parameters
fprintf('Standard');
tic;
    EM(map, d1n, 'stretch', 0.01, 'bend', 0.1);
toc
map1 = map;

% Create and initialise map
fprintf('creation');
map = rect2DMap(30, 30);
tic;
init(map, d1n, 'pci');
toc
fprintf('Trains');
tic;
    EM(map, d1n, 'stretch', 0.01, 'bend', 1);
toc
tic;
    EM(map, d1n, 'stretch', 0, 'bend', 1);
toc
tic;
    EM(map, d1n, 'stretch', 0, 'bend', 0.5);
toc
tic;
    EM(map, d1n, 'stretch', 0, 'bend', 0.1);
toc
map2 = map;

fprintf('Figures');
for k=1:2
    if k == 1
        map = map1;
    else
        map = map2;
    end
    tic;
    drawMapInt(map, d1n, 0, 'lineWidth', 0.5, 'nodeMarker', 'h');
    toc
end
for k=1:2
    if k == 1
        map = map1;
    else
        map = map2;
    end
    tic;
    drawMapInt(map, d1n, 1, 'lineWidth', 0.5, 'nodeMarker', 'h');
    toc
end
for k=1:2
    if k == 1
        map = map1;
    else
        map = map2;
    end
    tic;
    drawMapInt(map, d1n, 2, 'lineWidth', 0.5, 'nodeMarker', 'h');
    toc
end