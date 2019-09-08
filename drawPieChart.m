function drawPieChart(x, y, radius, Proportions, Colors)
%drawPieChart draw one pie chart with specified properties:
%   x and y are coordinates of centre of circle
%   radius is radius of circle
%   Proportions is vector of counts of cases for each color
%   Colors is vector of colors

    % Transform counts to proportions
    prop = Proportions / sum(Proportions);
    % Starts from zero angle
    angle = 0;
    % Draw chart
    if size(prop, 2) == 1
        % One class chart
        for i = 1:size(Proportions, 2)
            color = Colors(i,:);
            plot_circle(x, y, radius, color);
        end
    else
        % Several classes case
        for i=1:size(Proportions,2)
            th1 = angle;
            th2 = angle+prop(i)*2*pi;
            color = Colors(i,:);
            plot_arc(th1,th2,x,y,radius,color);
            angle = th2;
        end
    end
end

function plot_circle(x,y,r,color)
% Plot a circular arc as a pie wedge.
% a is start of arc in radians,
% b is end of arc in radians,
% (h,k) is the center of the circle.
% r is the radius.
% Try this:   plot_arc(pi/4,3*pi/4,9,-4,3)
% Author:  Matt Fig
    t = 0:0.05:2*pi;
    xp = r*cos(t) + x;
    yp = r*sin(t) + y;
    %x = [x h x(1)];
    %y = [y k y(1)];
    fill(xp,yp,color);
end

function plot_arc(a,b,h,k,r,color)
% Plot a circular arc as a pie wedge.
% a is start of arc in radians, 
% b is end of arc in radians, 
% (h,k) is the center of the circle.
% r is the radius.
% Try this:   plot_arc(pi/4,3*pi/4,9,-4,3)
% Author:  Matt Fig
    t = linspace(a,b);
    x = r*cos(t) + h;
    y = r*sin(t) + k;
    x = [x h x(1)];
    y = [y k y(1)];
    fill(x, y, color);
%     axis([h-r-1 h+r+1 k-r-1 k+r+1])
end