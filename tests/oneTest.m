function times = oneTest(data, sizes, stret, bend, optEM, optDraw)
%Provide one test for specified agguments
%
%Inputs:
%   data is n-by-m data matrix with n objects and m attributes.
%   sizes is two element vector with sizes of map. 
%   stret is k-element vector of stretch modulos.
%   bend is k-element vector of banding modulos.
%   optEM is cell arrays of arguments for EM.
%   optDraw is cellarray of arguments for drawMap and drawMapInt

    nSteps = length(stret);
    times = zeros(1 + nSteps, 1);

    %Create and initialise map and 
    map = rect2DMap(sizes(1), sizes(2));
    tic;
    init(map, data, 'pci');
    times(1) = toc;
    fprintf('Initialisation %g\n', times(1));
    
    % Draw initial map
    drawAll(map, data, 0, optDraw{:});
    
    for step = 1:nSteps
        % Training
        tic;
        EM(map, data, 'stretch', stret(step), 'bend', bend(step), optEM{:});
        times(step + 1) = toc;
        fprintf('Iteration %d %g\n', step, times(step + 1));
        drawAll(map, data, step, optDraw{:});
    end
end

function drawAll(map, data, iter, varargin)
    iter = num2str(iter);
    tit = ['Iteration ', iter];
    % Draw the map to see distribution
    drawMap(map, data, 'nodeMarker', 'none', 'lineWidth', 0.5,...
        varargin{:});
    title(tit);
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng', '-noui', '-loose', ['tmp\\', iter, '-1.png']);
    % Draw internal map with projection onto nodes
    drawMapInt(map, data, 0, 'lineWidth', 0.5, varargin{:});
    title(tit);
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng', '-noui', '-loose', ['tmp\\', iter, '-2.png']);
    % Draw internal map with projection onto edges
    drawMapInt(map, data, 1, 'nodeMarker', 'none', 'lineWidth', 0.5,...
        varargin{:});
    title(tit);
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng', '-noui', '-loose', ['tmp\\', iter, '-3.png']);
    % Draw internal map with projection onto edges
    drawMapInt(map, data, 1, 'nodeMarker', 'none', 'lineWidth', 0.5,...
        varargin{:});
    title(tit);
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng', '-noui', '-loose', ['tmp\\', iter, '-4.png']);
    close all;
end