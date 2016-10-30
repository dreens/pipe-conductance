%% Pipe Conductance
% Find the conductance of a long pipe with arbitrary cross section.
%
% Inputs:
% L - Length of the pipe, in centimeters.
% m - Average mass of a particle in the gas in kg. (4.65e-26 for N2)
% T - Temperature of the rarified gas in Kelvin. (293 K for room temp)
% pipeType - a string specifying the type of cross section. Must be one of:
%  circle, square, rectangle, triangle, polygon, annulus, wedge, boundary, image
%
% Additional Inputs:
%  Depending on the pipeType, different inputs will be expected directly
% after the pipeType input string:
% circle: diameter
% square: sidelength
% rectangle: [length, width] (in either order)
% triangle: sidelengths - a length 3 array of triangle sidelengths.
% polygon: vertices - an Nx2 array of vertices of the polygon. Cannot self-intersect.
% annulus: [radius #1, radius #2] (either order)
% wedge: [radius, angle]
% boundary: boundary - an arbitrary cross section with multiple disjoint components
%  or internal exclusions, specified as a 1D arraylist with Nx2 array elements, each
%  representing a connected component of the boundary as a list of ordered points
%  along the boundary, with the valid pipe area in which gas can flow to the right
%  of the points.
% image: {image, pixelsize} - image is a logical array representing a black and white
%  image of an arbitrary cross section, black is the inside of the pipe where gas 
%  can flow. pixelsize gives the size of a single pixel of the image in centimeters.
%
% Optional Inputs. After the above, property - value pairs can be specified:
% 'Velocity' - override the gas velocity calculated from temperature and gas mass 
%   with your own velocity specified in meters / second. m and T inputs must still
%   be provided but will not be used.
%
% 'Density' - specify how many points to break the perimeter into. Default is 500.
%   Be careful making this too large, as computation time is cubic in this parameter.
%
% 'Convex' - specify 'true' to force the program to assume a convex cross section.
%   Applies to polygon, boundary, and image only.
% 
% 'ViewImage' - specify 'true' to see a figure related to the computation while it
%   is performed. Will worsen performance, but can be useful for verifying expected
%   behavior for image, polygon, or boundary based cross sections.
%
% 'Units' - specify 'in', 'm', 'mm', 'ft', 'cm' (default) and pipe length and 
%   cross sectional parameters will be interpreted in those units.
%
% Outputs:
% c - Conductance of the pipe in the rarified gas regime (i.e. dilute
% molecular flow) measured in Liters / second.
function c = pipeConductance(L, m, T, pipetype, pipeargs, varargin)

% first lets interpret the inputs.
% We first find any string inputs in varargin and extract them with their
% value pairs.
properties = struct(...
    'Velocity',nan,...
    'Density',500,...
    'Convex',false,...
    'ViewImage',false,...
    'Units','cm');

assert(~mod(length(varagin),2),['Odd number of arguments after pipeargs.'...
    ' After pipeargs all further arguments must come as Name,Value pairs.'...
    ' Make sure pipeargs is a single input, wrapped as a list or'...
    ' arraylist as appropriate.'])

for i=1:2:length(varargin)
    name = varargin{1};
    assert(ischar(name),'Name argument in name value pair must be a string')
    value = varargin{2};
    if isfield(possiblePairs,name)
        properties.(name) = value;
    else
        warning('pipeCond:inputWarn',['There is no ''%s'' property for '...
            'Pipe Conductance calculations. Property and corresponding '...
            'value not set.'],name)
    end
end

assert(ischar(mode),'pipetype must be a string')
switch(pipetype)
    case 'circle'
        diameter = pipeargs;
        th = linspace(0,2*pi,properties.Denisty+1);
        th = th(1:end-1)';
        boundary = {diameter/2*[cos(th), -sin(th)]};
    case 'square'
        sidelength = pipeargs;
        a = round(properties.Density/4);
        row = linspace(0,sidelength,a+1)';
        xp = [zeros(a,1) ; row(1:end-1) ; ones(a,1)*sidelength ; flipud(row(2:end))];
        yp = [row(1:end-1) ; ones(a,1)*sidelength ; flipud(row(2:end)) ; zeros(a,1)];
        boundary = {[xp yp]};
    case 'rectangle'
        side1 = pipeargs(1);
        side2 = pipeargs(2);
        a = round(side1/(side1+side2)*properties.Density/4);
        b = round(side2/(side1+side2)*properties.Density/4);
        row = linspace(0,side1,a+1)';
        col = linspace(0,side2,b+1)';
        xp = [zeros(b,1) ; row(1:end-1) ; ones(b,1)*side1 ; row(end:-1:2)];
        yp = [col(1:end-1) ; ones(a,1)*side2 ; col(end:-1:2) ; zeros(a,1)];
        boundary = {[xp yp]};
    case 'triangle'
        side1 = pipeargs(1);
        side2 = pipeargs(2);
        side3 = pipeargs(3);
        sidesort = sort(pipeargs);
        assert(sidesort(end) < sum(sidesort(1:2)),'pipeCond:triangleErr',...
            'Sidelengths %f, %f, and %f do not make a triangle.',...
            side1,side2,side3);
        angle12 = acos((side1^2+side2^2-side3^2)/(2*side1*side2));
        pvx = side2*sin(angle12);
        pvy = side1-side2*cos(angle12);
        a = round(side1/(side1+side2+side3)*properties.Density/3);
        b = round(side2/(side1+side2+side3)*properties.Density/3);
        c = round(side3/(side1+side2+side3)*properties.Density/3);
        x1 = zeros(a,1);
        y1 = linspace(0,side1,a);
        x2 = linspace(0,pvx,b);
        y2 = linspace(side1,pvy,b);
        x3 = linspace(pvx,0,c);
        y3 = linspace(pvy,0,c);
        boundary = [x1(1:end-1) x2(1:end-1) x3(1:end-1) ; ...
                    y1(1:end-1) y2(1:end-1) y3(1:end-1) ];
        boundary = boundary';
    case 'polygon'
        vertices = pipeargs;
        assert(size(vertices,2)==2,'pipeCond:polygonErr',...
            'Vertex array must have size Nx2, not Nx%d',size(vertices,2));
        assert(size(vertices,1)>2,'pipeCond:polygonErr',...
            'Vertex array must have at least 3 vertices');
        sidels = sqrt(sum(vertices.^2,2));
        tl = sum(sidels);
        sideps = round(sidels/tl)*properties.Density);
        
        


% First we add points along the border of the polygon, then we pass it to
% getAreaInertia.
    






end