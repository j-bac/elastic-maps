function data = createTestCurve( dim, N, A, R, M)
%Generate N dim dimensional data points arround the arc of circle of radius
%R from angle -A to A with uniform noise with each coordinate magnitude M.
%dim is number of coordinates for each point
%N is number of points
%A is angle of arc
%R is radius of circle
%M is magnitude of noise

    % Generate Noise
    data = (rand(N, dim) - 0.5) * M;
    % Generate Angles
    ang = (rand(N, 1) - 0.5) * (2 * A);
    % Add basic circle to noise
    data(:, 1) = data(:, 1) + R * cos(ang);
    data(:, 2) = data(:, 2) + R * sin(ang);
end

