function EM(map, data, varargin)
%EM is function to fit map to data
%   Syntax
%   EM( map, data )
%   EM( __, Name, Value )
%   
%Inputs:
%   map is an object of MapGeometry class or subclass of this class.
%   data is n-by-m matrix which is the data set to fit. Each row contains
%       m coordinates of one data point.
%   Name can be one of the following:
%       'type' is one of the following strings:
%          'hard' is hard map with stretch = 1 and bend = 1
%           'medium' is more flexible map with stretch = 0.7 and bend = 0.7
%           'soft' is soft map with stretch = 0.5 and bend = 0.5
%           If 'type', 'stretch' and 'bending' are omitted then 'medium' is
%           used. 
%       'stretch' is a positive numeric value which is the value of
%           stretching modulo or a function with syntax
%               val = stretch( epoch )
%           where epoch is number of epoch (see epoch definition below) and
%           val is the nonnegative stretching modulo to use on specified
%           epoch. Epochs are numerated from 1.
%           Default value corresponds to type 'medium'
%       'bend' is a positive numeric value which is the value of the
%           bending modulo or a function with syntax 
%               val = bend( epoch )
%           where epoch is number of epoch (see epoch definition below) and
%           val is the bending modulo to use on specified epoch. Epochs are
%           numerated from 1. 
%           Default value corresponds to type 'medium'
%       'weights' is n-by-1 vector of weights for data points. Weights must
%           be nonnegative.
%       'intervals', intervals serves to specify user defined intervals.
%           intervals is row vector. The first element must be zero. By
%           default is created by usage of 'number_of_intervals' and
%           'intshrinkage'. Maximal value M is calculated as maximum among
%           data points of distance from data points to the nearest node of
%           map after initiation. Then this value is multiplied by
%           'intshrinkage'. All other borders are calculated as r(i) =
%           M*i^2/p^2, where p is number_of_intervals'. Ignored if
%           'potential' is not specified.
%       'Number_of_intervals' specifies the number of intervals to
%           automatic interval calculation. Default value is 5. Ignored if
%           'potential' is not specified. 
%       'intshrinkage' is fraction of maximal distance from data points to
%           original map which is used for intervals shrinkage (see
%           argument delta in defineIntervals). Default value is 1 (no
%           shrinkage). Ignored if 'potential' is not specified.
%       'potential' is majorant function for PQSQ. L2 distance without
%           shrinkage is used if 'potential' is not specified.
%
% One epoch is fitting of map with fixed values of stretching and bending
% modulo. This process can include several iterations of two step
% algorithm:
%   1. associate each data point with nearest node.
%   2. recalculate node position.
% Process of map fitting is stopped if new values of stretching and bending
% modulo are the same as on previous epoch OR if both stretching and
% bending modulo are zero.

    % Check the number of input attributes and types of the two first
    % attributes.
    if nargin < 2
        error('At least map and data must be specified');
    end    
    if ~isa(map,'MapGeometry')
        error('Incorrect type of the "map" argument, it must be MapGeometry');
    end
    if ~ismatrix(data) || ~isnumeric(data)
        error('Incorrect type of the "data" argument, data must be a matrix');
    end
    
    % Data preprocessing
    data = map.preprocessData(data);
    
    % Get sizes of data
    [n, dim] = size(data);

    % Default values of customisable variables
    strFun = @constStretch;
    constStretching = 0.7;
    bendFun = @constBend;
    constBending = 0.7;
    weights = [];
    func = [];
    intervals = [];
    nInt = 5;
    delta = 1;
    
    % Decode varargin
    for i=1:2:length(varargin)
        if strcmpi(varargin{i}, 'type')
            switch lower(varargin{i + 1})
                case 'hard'
                    strFun = @constStretch;
                    constStretching = 1;
                    bendFun = @constBend;
                    constBending = 1;
                case 'medium'
                    strFun = @constStretch;
                    constStretching = 0.7;
                    bendFun = @constBend;
                    constBending = 0.7;
                case 'soft'                    
                    strFun = @constStretch;
                    constStretching = 0.5;
                    bendFun = @constBend;
                    constBending = 0.5;
                otherwise
                    error('Incorrect value for type argument');
            end
        elseif strcmpi(varargin{i}, 'stretch')
            tmp = varargin{i + 1};
            if isa(tmp, 'function_handle')
                strFun = tmp;
            else
                strFun = @constStretch;
                constStretching = tmp;
            end
        elseif strcmpi(varargin{i}, 'bend')
            tmp = varargin{i + 1};
            if isa(tmp, 'function_handle')
                bendFun = tmp;
            else
                bendFun = @constBend;
                constBending = tmp;
            end
        elseif strcmpi(varargin{i}, 'weights')
            weights = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'intervals')
            intervals = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'Number_of_intervals')
            nInt = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'intshrinkage')
            delta = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'potential')
            func = varargin{i + 1};
        else
            if ischar(varargin{i})
                error(['Wrong name of argument "', varargin{i}, '"']);
            else
                error(['Wrong name of argument "', num2str(varargin{i}), '"']);
            end
        end
    end

    
    % Check type and length of weights
    if isempty(weights)
        weights = ones(n, 1);
    else
        %Weights must be a vector of nonnegative finite reals with at least two
        %values greater than zero and with number of elements equal to number
        %of rows in X. 
        if ~isreal(weights) || ~isfinite(weights) || sum(weights<0)>0 ||...
                sum(weights>0)<2 || size(weights, 1) ~= n ||...
                size(weights, 2) ~= 1
            error(['Incorrect value for argument "Weights". It must be ',...
                'a column vector of nonnegative finite reals with at',...
                'least two values greater than zero and with number',...
                ' of elements equal to number of rows in data.']);
        end
        weights = weights(:);
    end

    % Define total weights
    TotalWeight = sum(weights);
    weigh = weights;
    pFunc = [];
    
    % Analyse PQSQ request
    if ~isempty(func)
        %Func must be function handler
        if ~isa(func,'function_handle')
            error(['Incorrect value in "potential" argument.'...
                ' It must be function handler']);
        end

        if isempty(intervals)
            %Function has to create intervals by automatic way
            %nInt must be positive integer scalar
            if ~isreal(nInt) || ~isfinite(nInt) || nInt < 1
                error(['Incorrect value of "number_of_intervals" argument' ...
                    'It must be positive integer scalar']);
            else
                nInt = floor(nInt);
            end
            %delta has to be positive real scalar
            if ~isreal(delta) || ~isfinite(delta) || delta < 0
                error(['Incorrect value of "intshrinkage" argument' ...
                    'It must be positive real scalar']);
            end
            
            pFunc = definePotentialFunction(map.getDisp(), nInt, func, delta);
        else
            %intervals must contains non negative values in ascending order.
            %The first value must be zero.
            if intervals(1)~=0 || ~all(isfinite(intervals)) ...
                    || any((intervals(2:end)-intervals(1:end-1))<=0)
                error(['Incorrect values in argument intervals: intervals must'...
                    ' contains finite non negative values in ascending order.'...
                    ' The first value must be zero.']);
            end
            pFunc.intervals = [intervals(:)', Inf(1)];
            [pFunc.A, pFunc.B] = ...
                computeABcoefficients(intervals, func);
        end
    end
    
    %Get initial state of nodes
    nodes = map.getMappedCoordinates();
    if size(nodes, 2) ~= dim
        error('Dimensions of mapped nodes and data must be the same');
    end
    N = size(nodes, 1);

    %Form matrices B and C
    tmp = map.getLinks();    
    B = diag(accumarray(tmp(:), 1, [N, 1]));
    tmp = accumarray(tmp, 1, [N, N]);
    B = B - tmp - tmp';
    
    tmp = map.getRibs();
    C = diag(accumarray([tmp(:, 1); tmp(:, 3)], 1, [N, 1])...
        +accumarray(tmp(:, 2), 4, [N, 1]));
    w = accumarray(tmp(:, [1, 3]), 1,[N, N]);
    tmp = accumarray([tmp(:, 1:2); tmp(:, 2:3)], 2, [N, N]);
    C = C + w + w' - tmp - tmp';

    % Start iterative process
    epoch = 1; % Number of iteration
    ass = zeros(n, 1); % Initial associations. It is impossible combination
    qInd = zeros(n, 1);
    % Get initial modulo
    stretch = strFun(epoch);
    bend = bendFun(epoch);
    while true
        % Save old associations and q indices.
        oldAss = ass;
        oldQInd = qInd;
        % Find new associations
        [dist, ass] = associate(map, nodes, data);
        % Find indeces for PQSQ if required
        if ~isempty(pFunc)
            [~,qInd]=histc(dist,pFunc.sqint);
            weigh = weights .* pFunc.A(qInd)';
        end

        % If nothing is changed then we have end of epoch
        if all(oldAss == ass) && all(oldQInd == qInd)
            epoch = epoch + 1;
            tmp = strFun(epoch);
            tmp1 = bendFun(epoch);
            if tmp == 0 && tmp1 == 0
                break;
            end

            if abs(tmp - stretch) + abs(tmp1 - bend) == 0
                break;
            end
            stretch = tmp;
            bend = tmp1;
        end
        
        % Form matrix A
        % For further robust and so on we consider possibility of zeros in
        % ass and create dummy element
        ass = ass + 1;
        % Calculate number of points for each node
        tmp = accumarray(ass, weigh, [N + 1, 1]) ;
        % Normalise and remove dummy element
        NodeClusterRelativeSize = tmp(2:end) / TotalWeight;
        % Create centroids
        NodeClusterCenters = zeros(N + 1, dim);
        for k = 1:dim
            NodeClusterCenters(:, k) =accumarray(ass, data(:, k) .* weigh, [N + 1, 1]) / TotalWeight
        end
        % Remove dummy element
        NodeClusterCenters = NodeClusterCenters(2:end,:);
        
        % form SLAE
        SLAUMatrix = diag(NodeClusterRelativeSize) + stretch * B + bend * C;
        nodes = SLAUMatrix \ NodeClusterCenters;
    
        % Restore ass
        ass = ass - 1;
    end
    
    % Put new nodes into map
    map.putMapped(nodes);
    
    function stretch = constStretch( ~ )
        stretch = constStretching;
    end

    function bend = constBend( ~ )
        bend = constBending;
    end
end

function potentialFunction = definePotentialFunction( x,...
    number_of_intervals, potential_function_handle, delta )
%definePotentialFunction defines "uniform in square" intervals for trimming
%threshold x and specified number_of_intervals.
%   x is upper boundary of the interval last but one.
%   number_of_intervals is required number of intervals.
%   potential_function_handle is function handler for coefficients
%       calculation.
%   delta is coefficient of shrinkage which is greater than 0 ang not
%       greater than 1.
%Output argument potentialFunction is structure with three fields:
%   intervals is matrix m-by-number_of_intervals. Each row contains
%       number_of_intervals values of thresholds for intervals and one
%       additional value Inf
%   A and B are the m-by-number_of_intervals matrices with quadratic
%       functions coefficients

    if nargin < 4 
        delta = 1;
    end
    
    p = number_of_intervals - 1;
    
    %intervals is the product of row and maximal coefficient multiplied by delta:
    intervals = (x * delta) * ((0:p) / p) .^ 2;
    potentialFunction.intervals = [intervals, Inf(1)];
    potentialFunction.sqint = potentialFunction.intervals .^ 2;
    [potentialFunction.A,potentialFunction.B] = ...
        computeABcoefficients(intervals, potential_function_handle);
end

function [A,B] = computeABcoefficients(intervals, potential_function_handle)
%PQSQR_computeABcoefficients calculates the coefficients a and b for
%quadratic fragments of potential function.
%   intervals is the 1-by-K matrix of intervals' boundaries without final
%       infinit boundary.
%   potential_function_handle is a handle of majorant function.

    %Get dimensions of intervals
    p = size(intervals,2);

    %Preallocate memory
    A = zeros(1,p);
    B = zeros(1,p);

    %Calculate value of function all boundaries
    pxk = potential_function_handle(intervals);
    sxk = intervals.^2;

    A(1:p-1) = (pxk(1:p-1)-pxk(2:p))./(sxk(1:p-1)-sxk(2:p));
    B(1:p-1) = (pxk(2:p).*sxk(1:p-1)-pxk(1:p-1).*sxk(2:p))./...
        (sxk(1:p-1)-sxk(2:p));
    B(p) = pxk(p);
end
