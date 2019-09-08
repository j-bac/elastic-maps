function data = createTestSurf( dim, N, A, R, M)
%Generate N dim dimensional data points arround the fragment of shphere of
%radius R with two polar angles: azimutal angle from -A to A and polar
%angle from pi to pi. Uniform noise with each coordinate magnitude M. 
%dim is number of coordinates for each point
%N is number of points
%A is azimutal angle
%R is radius of sphere
%M is magnitude of noise

    % Generate Noise
    data = (rand(N, dim) - 0.5) * M;
    % Generate Angles
    bng = (rand(N, 1) - 0.5) * (2 * A);
    ang = (rand(N, 1) - 0.5) * pi;
    % Add basic circle to noise
    data(:, 1) = data(:, 1) + R * cos(ang) .* sin(bng);
    data(:, 2) = data(:, 2) + R * sin(ang) .* sin(bng);
    data(:, 3) = data(:, 3) + R * cos(bng);
end

