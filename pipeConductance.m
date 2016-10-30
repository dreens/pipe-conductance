%% Pipe Conductance
% Find the conductance of a long pipe with arbitrary cross section.
%
% Inputs:
% Peri - perimeter of polygon as an Nx2 list of vertices, in centimeters.
%        Disjoint polygons are also possible- 
% PipeL - Length of the pipe, also in centimeters.
% GasMass - Mass of a particle making up the gas in kg. (4.65e-26 for N2)
% Temp - Temperature of the rarified gas in Kelvin. (293 K for room temp)
% 
% Optional Inputs:
% density - specify how many segments to break the perimeter into. Default
% is 1000. Runtime increases as the cube of this parameter, so be careful
% specifying a larger number.
%
% Outputs:
% Cond - Conductance of the pipe in the rarified gas regime (i.e. dilute
% molecular flow) measured in Liters / second.
function c = pipeConductance(Peri, PipeL, GasMass, Temp, varargin)

% First we add points along the border of the polygon, then we pass it to
% getAreaInertia.
    






end