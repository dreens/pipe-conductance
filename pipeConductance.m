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
% rectangle: length, width (in either order)
% triangle: sidelengths - a length 3 array of triangle sidelengths.
% polygon: vertices - an Nx2 array of vertices of the polygon. Cannot self-intersect.
% annulus: radius #1, radius #2 (either order)
% wedge: radius, angle
% boundary: boundary - an arbitrary cross section with multiple disjoint components
%  or internal exclusions, specified as a 1D arraylist with Nx2 array elements, each
%  representing a connected component of the boundary as a list of ordered points
%  along the boundary, with the valid pipe area in which gas can flow to the right
%  of the points.
% image: image, pixelsize - image is a logical array representing a black and white
%  image of an arbitrary cross section, black is the inside of the pipe where gas 
%  can flow. pixelsize gives the size of a single pixel of the image in centimeters.
%
% Optional Inputs. After the above, name - value pairs can be specified:
% 'Velocity' - override the gas velocity calculated from temperature and gas mass 
%   with your own velocity specified in meters / second. m and T inputs must still
%   be provided but will not be used.
%
% 'Density' - specify how many points to break the perimeter into. Default is 1000.
%   Be careful making this too large, as computation time is cubic in this parameter.
%
% 'Convex' - specify 'true' to force the program to assume a convex cross section.
%   Applies to polygon, boundary, and image only.
% 
% 'ViewImage' - specify 'true' to see a figure related to the computation while it
%   is performed. Will worsen performance, but can be useful for verifying expected
%   behavior for image, polygon, or boundary based cross sections.
%
% 'Length' - specify 'in', 'm', 'mm', 'ft', 'cm' (default) and pipe length and 
%   cross sectional parameters will be interpreted in those units.
%
% Outputs:
% c - Conductance of the pipe in the rarified gas regime (i.e. dilute
% molecular flow) measured in Liters / second.
function c = pipeConductance(L, m, T, pipetype, varargin)

% first lets interpret the inputs.
assert(ischar(mode),'pipetype must be a string')
switch(pipetype)



% First we add points along the border of the polygon, then we pass it to
% getAreaInertia.
    






end
