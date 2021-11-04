function Position3D = RandomPosition3D(Dimension,Radius,Density)

 
% Dimension = [4096,4096,300];
% % Height = 500;
% Radius = 75;
% Density = 0.7;
 
% INPUT:
%   nWant:  Number of spheres
%   Dimension:  Dimension of 3d box as [1 x 3] double vector
%   Radius: Radius of spheres
% OUTPUT:
%   P:      [nWant x 3] matrix, centers

nWant  = 1;
P      = zeros(nWant, 3);
R2     = (2 * Radius) ^ 2;   % Squared once instead of SQRT each time
W      = Dimension - [2*Radius, 2*Radius,2*Radius]; % Avoid interesction with borders
iLoop  = 1;                  % Security break to avoid infinite loop
nValid = 0;
while iLoop < 1e6
  newP = rand(1, 3) .* W;
  % Auto-expanding, need Matlab >= R2016b. For earlier versions:
  % Dist2 = sum(bsxfun(@minus, P(1:nValid, :), newP) .^ 2, 2);
  Dist2 = sum((P(1:nValid, :) - newP) .^ 2, 2);
  if all(Dist2 > R2)
    % Success: The new point does not touch existing sheres:
    nValid       = nValid + 1;  % Append this point
    P(nValid, :) = newP;
  end
  iLoop = iLoop + 1;
end
% Stop if too few values have been found:
if nValid < nWant
  error('Cannot find wanted number of points in %d iterations.', iLoop)
end
[column, row] = size(P);
MaxNumPart = column;
ss = randperm(MaxNumPart);
Prandom = P(ss(1:floor(Density*MaxNumPart)),:);
Prandom(:,1) = Prandom(:,1) - (Dimension(1)-2*Radius)*1/2;
Prandom(:,2) = Prandom(:,2) - (Dimension(2)-2*Radius)*1/2;
Prandom(:,3) = Prandom(:,3) - (Dimension(3)-2*Radius)*1/2;
Position3D = Prandom;

% % Plot the spheres in 3D space
% figure 
% axes('NextPlot', 'add', ...
%     'XLim', [-0.5*Dimension(1), 0.5*Dimension(1)], 'YLim', [-0.5*Dimension(2), 0.5*Dimension(2)], 'ZLim', [-0.5*Dimension(3), 0.5*Dimension(3)]);
% view(3);
% [X, Y, Z] = sphere();
% for k = 1:(Density*MaxNumPart)
%     surf(X * Radius + Prandom(k, 1), Y * Radius + Prandom(k, 2), Z * Radius + Prandom(k, 3));
% end
end
