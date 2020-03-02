function [circData, circHandles] = Randomposition(N, mindist)
% Create a figure 2D circles (or 'bubbles') with random centers and varying radii.
% Paramters are set and explained at the top of the code.
% Ouputs include
%   * circData: [m x n] matrix of m circles. Columns 1 & 2 are the [x,y] centers
%       of each cirlce. Column 3 is the radius for each circle.
%   * circHandles: [m x 1] vector of handles for each line object / circle.
% Examples of how to use outputs
%   * Calculate the area for each circle: circArea = pi * circData(:,3).^2;
%   * Change color and line width of circle edges:  set(circHandles, 'color', 'r', 'LineWidth', 3)
%   * Fill circles with color: %
%         c = jet(size(circData,1)); %choose colors
%         for i = 1:size(circData,1)
%             rectangle('Position', [circData(i,1:2)-circData(i,3), circData(i,[3,3])*2], 'Curvature',[1,1], 'FaceColor',c(i,:))
%         end
%   * Count number of circles for each circle size:
%       [g, radius] = findgroups(circData(:,3));
%       count = splitapply(@length, circData(:,3), g);
%       table(radius, count)
% Details (other than what's explained in the parameter section)
%   * For each circle size (defined by user), the number of circles that are
%       attempted is a proportion between the circle's area and the frame's area.
%       For 'n' circle sizes, the total area of the frame is evenly split into 'n'
%       portion and the number of circles per circle size is merely framePortion / circleArea.
%   * Circles are drawn in order of circle size from largest to smallest. The distance
%       between circles is calculated by the distance between all circle centers minus
%       the radii.  That gives the distance between all circle edges.
%   * Since circle centers are chosen randomly, we check that the edges won't overlap
%       with existing circles.  If they do, we choose another set of random centers
%       and that process repeats in a while loop until either all circle centers
%       do not overlap or we meet a maximum iteration (set by user).  A warning
%       message appears if the while loop expired without drawing all expected
%       circles.
%
% This function was written for a Matlab forum question:
% https://www.mathworks.com/matlabcentral/answers/446114-non-overlapping-random-circles
% 02/21/2019
% Adam Danz

%% set parameters
% Frame size (arbitrary units): size of axes, centered at (0,0).
fs = [floor(N-1*mindist), floor(N-1*mindist)]; %[width, height]

% Circle sizes: Select the minimum and maximum radii and
%   the number of unique radii.
%   Radii are linearly spaced. ie: linspace(cs(1), cs(2), nSizes).
cs = [floor(mindist*0.5), floor(mindist*0.5)]; %[smallest, largest] radius
nSizes = 1;    %number of circle sizes

% max iterations: how many attempts should be made to find circle
%   locations that do not overlap with other circles?
maxIt = 100;

% Decide whether circles shoule be entirely inside of the frame
% 0 = Cirlces edges can expand outside of the frame
% 1 = Cirlces must be entirely inside of the frame
% 2 = Circle edges that extend beyond the frame are cut off
edgeType = 1;

% Density of circles:
% a value greater than 0 but less than or equal to 1.
% As the value approaches 0, less circles will be drawn.
density = 1;

% Allow overlap (true or false)
% When true, circles may overlap.
allowOverlap = false;

%% Add soap
% Create figure
% figure();
% axes();
% hold on
% axis equal %set aspect ratio to 1:1
% xlim(fs(1)/2 * [-1.05, 1.05]) %with some extra space
% ylim(fs(2)/2 * [-1.05, 1.05]) %with some extra space

% shrink frame so max cirlce is always inside
switch edgeType
    case 0  % Cirlce edges can expand outside of the frame
        fs2 = fs;
        clearBorder = false;
    case 1  % Cirlce must be entirely inside of the frame
        fs2 = fs - max(cs)*2;
        clearBorder = false;
    case 2  % Circle edges end at the frame
        fs2 = fs;
        clearBorder = true;
    otherwise
        error('Did not recognize edgeType.')
end

% determine minimum distance between circles allowed
if allowOverlap
    minDist = -inf;
else
    minDist = 0;
end

% determine circle sizes (largest to smallest)
r = linspace(cs(2), cs(1), nSizes);

% determine approximate number of circles to draw for each size
frameArea = fs(1) * fs(2);
circAreas = pi * r.^2;
d = ceil(ceil((frameArea / nSizes) ./ circAreas) * density); %Old: (2:nSizes+1).^2;

% Throw error if largest circle is larger than frame
if max(circAreas) > frameArea
    error('The area of the largest circle (%.1f) is larger than the frame''s area (%.1f)', max(circAreas), frameArea)
end

% Loop through each circle size
circData = []; %[xCenter, yCenter, radius] of each drawn circle
h = cell(nSizes,1); % handles to the line objects for each circle
% wb = waitbar(0,'initializing...','name',mfilename); % waitbar
for i = 1:nSizes
    % Reset circle & iteration count
    cCount = 0;
    iCount = 0;
    xRand = zeros(d(i),1);
    yRand = xRand;
    isOK = false(size(xRand));
    % Keep drawing random coordinates until either all circs
    % are drawn or we reach max number of attempts
    while cCount < d(i) && iCount < maxIt
        % Randomly choose center coordinates
        xRand(~isOK) = (rand(d(i)-cCount,1) - 0.5) * fs2(1);
        yRand(~isOK) = (rand(d(i)-cCount,1) - 0.5) * fs2(2);
        % determine if new circles overlap with others or themselves
        xyr = [circData; [xRand, yRand, repmat(r(i), size(yRand))]];
        radMat = repmat(xyr(:,3),1,size(xyr,1)) + repmat(xyr(:,3).',size(xyr,1),1); %changed 190901 to work with r2016a
        dist = tril(squareform(pdist(xyr(:,1:2))) - radMat, -1);
        if isempty(dist)
            isOK = true; %when xyr only has 1 row
        else
            isOK = all(dist(size(circData,1)+1:end,:) >= minDist, 2);
        end
        cCount = sum(isOK);  %cirlce count for current radius
        iCount = iCount + 1; %iteration count
%         % Update waitbar     %waitbar
%         if ishghandle(wb)
%             waitbar(max(iCount/maxIt,cCount >= d(i)),wb,sprintf('Trying to find space for circles with radius = %.2f (min:%.2f)',r(i),r(end)));
%         end
     end
    % If we had to quit searching, throw warning.
    if iCount >= maxIt
        warning('Max iteration reached.  Only %d of the requested %d circles drawn for radius %.3f', cCount, d(i), r(i))
    end
    % Store all final circle data
    circData = [circData; [xRand(isOK), yRand(isOK), repmat(r(i), sum(isOK), 1)]];
%     % Draw circles
%     if any(isOK)
%         h{i} = drawcircles([xRand(isOK), yRand(isOK), repmat(r(i), sum(isOK))], clearBorder, fs);
%     end
    
end
% % Remove waitbar
% if ishghandle(wb)
%     delete(wb)
% end
% 
% % Draw frame
% rectangle('position', [-fs/2, fs], 'LineWidth', 2)
% circHandles = [h{:}]';

% function h = drawcircles(xyr, clearBorder, fs)
% % Draw circle given center and radius
% ang=0:0.01:2*pi;
% xp=xyr(:,3)*cos(ang) + repmat(xyr(:,1),1,numel(ang)); %changed 190901 to work with r2016a
% yp=xyr(:,3)*sin(ang) + repmat(xyr(:,2),1,numel(ang)); %changed 190901 to work with r2016a
% if clearBorder
%     % remove data outside of frame, if requested
%     xp(abs(xp) > fs(1)/2) = NaN;
%     yp(abs(yp) > fs(2)/2) = NaN;
% end
% h = plot(xp', yp', 'k')';
