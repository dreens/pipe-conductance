%% getPerimeter description. 
% A bw image of a masked area is converted to a list of perimeter points
% that can be passed to getAreaInertia. This involves scaling all of the
% perimeter points from unitless integers to floats representing position
% in a 2D plane measured in centimeters. This also involves assuring the
% handedness of each component of the boundary, so that the area is always
% on the right. The bulk of the work is done by MATLAB's bwboundaries.
%
% sp should give the width in centimeters of a single pixel. Care should be
% taken to provide crossArea in a square manner. Obviously it would be
% possible to implement a different centimeter/pixel ratio for each
% direction, but this is more error-prone.
function B = getPerimeter( sp , crossArea )

[B,~,N,~] = bwboundaries(crossArea);

for k=1:length(B)
    
    % first scale all the points
    a = B{k}*sp;

    % bwboundaries returns the outer boundaries first, then the inner. N
    % gives the number of outer boundaries only, so only inner boundaries
    % pass this conditional. bwboundaries gives all boundaries clockwise,
    % so outer already satisfy "area on the right" but inner need
    % reversing.
    if k>N
        a = flipud(a);
    end
    
    % Here we throw away some of the points. bwboundaries forces points
    % onto a grid, throwing some away is a rudimentary form of smoothing.
    % We could do better by applying a smoothing spline, but not now.
    % Throwing away points is also important for run-time, since
    % getAreaInertia is cubic in boundary length.
    a = a(1:50:end-1,:);
    
    B{k} = a;

end

end