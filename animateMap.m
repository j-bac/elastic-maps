function animateMap(fig, filePrefix, startAngles, steps)
%animateMap rotates figure fig from the inuitial position startAngles with
%increment steps until the full rotation.
%Inputs:
%   fig is handle of axes to rotate. For example fig = cga; or 
%       fig = subplot(...);.
%   filePrefix is text constant for file names. In each iteration k file
%       with name [filePrefix, sprintf(%3d, k), '.png'] is created.
%   startAngles is one or two element vector. If startAngles contains one
%       element then it is Azimuth angle, and elevation angle is 30.
%       Otherwise the first element is Azimuth vector and the second
%       element is Elevation vector. For detailed angle description see
%       function view in standard Camera Views documentation. 
%       Default value is [-37.5, 30]. 
%   steps is one or two dimensional vector with angles incrementation. If
%       the second element is omitted the it is equal 0. Default value is
%       [10, 0];
%
    % Sanity check of parameters
    if nargin < 2
        error('At least fig and filePrefix must be presented');
    end
    if ~isgraphics(fig) || ~strcmp(get(fig,'type'),'axes')
        error('Wrong type of the first argument. It must be axes handle');
    end
    if nargin < 3 || isempty(startAngles)
        startAngles = [-37.5, 30];
    end
    if numel(startAngles) == 1 
        startAngles = [startAngles, 30];
    end
    startAngles = startAngles(:)';
    if startAngles(2) > 180 
        startAngles = startAngles - 360;
    end
    if nargin < 4 || isempty(steps)
        steps = [10, 0];
    end
    if numel(steps) == 1 
        steps = [steps, 0];
    end
    steps = steps(:)';
    
    % Activate specified graphics object
    f = get (fig, 'Parent');
    set(0, 'currentfigure', f);
    set(f, 'currentaxes', fig);
    
    % Fix sizes
    set(gcf, 'PaperPositionMode', 'auto');
    axis vis3d;

    % Calculate number of steps
    if steps(1) ~= 0
        K = 360 / abs(steps(1));
    else
        K = 0;
    end
    if steps(2) ~= 0
        KK = 360 / abs(steps(2));
    else
        KK = 0;
    end
    if K < KK
        K = KK;
    end
    if K == 0 
        error('Wrong steps parameter');
    end
    % Now we ready to start rotation
    view(fig, angleNorm(startAngles));
    k = 0;
    print(gcf,'-dpng', '-noui', '-loose',...
        [filePrefix, sprintf('%03d', k), '.png']);
    for k = 1:K
        startAngles = startAngles + steps;
        view(angleNorm(startAngles));
        print(gcf,'-dpng', '-noui', '-loose',...
            [filePrefix, sprintf('%03d', k), '.png']);
    end
end

function angl = angleNorm(angles)
% Normalises angles for arbitrary Elevation angle.
    angl = angles;
    if angles(2) > 90 && angles(2) < 270
        angl(1) = angl(1) - 180;
    end
end