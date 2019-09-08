% Colours to draw are defined now
colours = ['b'; 'g'];

% Download database
load('breastCancer.mat');

if size(d1n, 1) > 286
    %We need to transpose matrix
    d1n = d1n';
end

% Create map wit 30 rows and 30 columns and initialise along Principal
% components
map = rect2DMap(30, 30);
init(map, d1n, 'pci');

% Draw the firs map to see distribution (map hidden by 'lineWidth', 0)
drawMap(map, d1n, 'classes', col, 'markColour', colours, 'lineWidth', 0);
saveFigures('figures/Distribution.png');

% Map fitting for standard parameters
EM(map, d1n, 'stretch', 0, 'bend', 0.1);

% Draw the density by density function
drawMapInt(map, d1n, 2, 'nodeMarker', 'none', 'lineWidth', 0.5,...
            'classes', col, 'markColour', colours, ...
            'coloring', 'density', 'flatColoring', true, 'ColorMap', parula);
saveFigures('figures/flatDensityC.png');
drawMapInt(map, d1n, 2, 'nodeMarker', 'none', 'lineWidth', 0.5,...
            'classes', col, 'markColour', colours, ...
            'coloring', 'density', 'flatColoring', false, 'ColorMap', parula);
saveFigures('figures/3DDensityC.png');
% Draw density by vector of ones
col1 = col;
col1(col1 > 0) = 1;
drawMapInt(map, d1n, 2, 'nodeMarker', 'none', 'lineWidth', 0.5,...
            'classes', col, 'markColour', colours, ...
            'coloring', col1, 'flatColoring', true, 'ColorMap', parula);
saveFigures('figures/flatDensityC2.png');
drawMapInt(map, d1n, 2, 'nodeMarker', 'none', 'lineWidth', 0.5,...
            'classes', col, 'markColour', colours, ...
            'coloring', col1, 'flatColoring', false, 'ColorMap', parula);
saveFigures('figures/3DDensityC2.png');
% Draw pseudo risk map
col1 = col;
col1(col1 == 2) = -1;
col1 = -col1;
drawMapInt(map, d1n, 2, 'nodeMarker', 'none', 'lineWidth', 0.5,...
            'classes', col, 'markColour', colours, ...
            'coloring', col1, 'flatColoring', true, 'ColorMap', parula);
saveFigures('figures/flatRiskC.png');
drawMapInt(map, d1n, 2, 'nodeMarker', 'none', 'lineWidth', 0.5,...
            'classes', col, 'markColour', colours, ...
            'coloring', col1, 'flatColoring', false, 'ColorMap', parula);
saveFigures('figures/3DRiskC.png');
