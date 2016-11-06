function c = pipeConductance(L, pipetype, pipeargs, varargin)
%PIPECONDUCTANCE Find the conductance of a pipe of arbitrary cross section.
%
%c = pipeConductance(L, pipetype, pipeargs) finds the conductance of a pipe
%   with length L in centimeters, with a cross section of type pipetype, a
%   string specifying the type of cross section, and with parameters
%   varying by pipetype but specified in the pipeargs input.
%
%   pipetype must be one of: circle, square, rectangle, triangle, polygon,
%   annulus, wedge, boundary, image
%
%   In the following, the pipeargs format corresponding to a given pipetype
%   is specified:
%   pipetype -> pipeargs         further description
%    circle      diameter
%
%    square      sidelength
%
%    rectangle   [length, width]  either order
%
%    triangle	 [a, b, c]        length 3 array of triangle's sidelengths
%
%    polygon     vertices         Nx2 array of vertices of the polygon
%
%    annulus     [r1, r2]         either order
%
%    wedge       [r, angle]
%
%    boundary    boundary         an arbitrary cross section with multiple
%                                 disjoint components or internal
%                                 exclusions, specified as a 1D arraylist
%                                 whose elements are Nx2 arrays, each
%                                 representing a connected component of the
%                                 boundary as a list of ordered points
%                                 along it, with the interior to the right.
%
%    image       {image, pixsize} image is a black and white image of an
%                                 arbitrary cross section, black is the
%                                 pipe interior. pixsize gives the size of
%                                 a single pixel of the image in
%                                 centimeters.
%
% c = pipeConductance(L, pipetype, pipeargs, 'Property', Value) specifies
%   additional arguments as Property, Value pairs. Valid properties are
%   listed below:
%
% 'Mass'    Specify the mass of the gas in AMU. Default is 29, which is the
%           average mass per molecule of "air" by which we mean three parts
%           molecular nitrogen and one part molecular oxygen.
%
% 'Temp'    Gas temperature in Kelvins. Default is 293 (20C, 68F).
%
% 'Gas'     Provide the type of gas as a string. This is just a shortcut
%           for specifying the mass. Valid types: H2, He, Ne, N2, Air, O2,
%           Ar, Kr, Xe, SF6, and many others, all expressed as chemical
%           formulas with capitalization as appropriate from the periodic
%           table.
%
% 'Velocity' Override the gas velocity calculated from temperature and gas
%           mass with your own velocity specified in meters / second. By
%           default, velocity is computed from 'Mass' and 'Temperature'.
%
% 'Density' Specify how many points to break the perimeter into. Default is
%           500. Be careful making this too large, as computation time is
%           cubic in this parameter for non-convex areas. 
%
% 'ViewImage' Default is false. If set true, a figure related to the
%           computation is shown. Will worsen performance, but can
%           be useful for verifying expected behavior for image, polygon,
%           or boundary based cross sections.
%
% 'Units' - Select 'in', 'm', 'mm', 'ft', 'cm' (default) and pipe length
%           and cross sectional parameters will be interpreted in those
%           units. Some other unusual units also implemented.
%
%
%
%
% Outputs:
% c -       Conductance of the pipe in the rarified gas regime (i.e. dilute
%           molecular flow) measured in Liters / second.

%% Interpret the inputs.
% We first find any string inputs in varargin and extract them with their
% value pairs.
props = struct(...
    'Mass',29,...
    'Temp',293,...
    'Gas',nan,...
    'Velocity',nan,...
    'Density',500,...
    'ViewImage',false,...
    'Units','cm');

assert(~mod(length(varargin),2),['Odd number of arguments after pipeargs.'...
    ' After pipeargs all further arguments must come as Name,Value pairs.'...
    ' Make sure pipeargs is a single input, wrapped as a list or'...
    ' arraylist as appropriate.'])

for i=1:2:length(varargin)
    name = varargin{i};
    assert(ischar(name),'Name argument in name value pair must be a string')
    value = varargin{i+1};
    if isfield(props,name)
        props.(name) = value;
    else
        warning('pipeCond:inputWarn',['There is no ''%s'' property for '...
            'Pipe Conductance calculations. Property and corresponding '...
            'value not set.'],name)
    end
end

assert(ischar(pipetype),'pipetype must be a string')
switch(pipetype)
    case 'circle'
        diameter = pipeargs;
        th = linspace(0,2*pi,props.Density);
        th = th';
        boundary = {diameter/2*[cos(th), -sin(th)]};
    case 'square'
        sidelength = pipeargs;
        a = round(props.Density/4);
        row = linspace(0,sidelength,a)';
        xp = [zeros(a,1) ; row ; ones(a,1)*sidelength ; flipud(row)];
        yp = [row ; ones(a,1)*sidelength ; flipud(row) ; zeros(a,1)];
        boundary = {[xp yp]};
    case 'rectangle'
        side1 = pipeargs(1);
        side2 = pipeargs(2);
        a = round(side1/(side1+side2)*props.Density/2);
        b = round(side2/(side1+side2)*props.Density/2);
        row = linspace(0,side1,a)';
        col = linspace(0,side2,b)';
        xp = [zeros(b,1) ; row ; ones(b,1)*side1 ; row(end:-1:1)];
        yp = [col ; ones(a,1)*side2 ; col(end:-1:1) ; zeros(a,1)];
        boundary = {[xp yp]};
    case 'triangle'
        side1 = pipeargs(1);
        side2 = pipeargs(2);
        side3 = pipeargs(3);
        sidesort = sort(pipeargs);
        assert(sidesort(end) < sum(sidesort(1:2)),'pipeCond:triangleErr',...
            'Sidelengths %g, %g, and %g do not make a triangle.',...
            side1,side2,side3);
        angle12 = acos((side1^2+side2^2-side3^2)/(2*side1*side2));
        pvx = side2*sin(angle12);
        pvy = side1-side2*cos(angle12);
        a = round(side1/(side1+side2+side3)*props.Density);
        b = round(side2/(side1+side2+side3)*props.Density);
        c = round(side3/(side1+side2+side3)*props.Density);
        x1 = zeros(1,a);
        y1 = linspace(0,side1,a);
        x2 = linspace(0,pvx,b);
        y2 = linspace(side1,pvy,b);
        x3 = linspace(pvx,0,c);
        y3 = linspace(pvy,0,c);
        boundary = [x1 x2 x3 ; ...
                    y1 y2 y3 ];
        boundary = {boundary'};
    case 'polygon'
        vertices = pipeargs;
        assert(size(vertices,2)==2,'pipeCond:polygonErr',...
            'Vertex array must have size Nx2, not Nx%d',size(vertices,2));
        assert(size(vertices,1)>2,'pipeCond:polygonErr',...
            'Vertex array must have at least 3 vertices');
        vertices = [vertices ; vertices(1,:)];
        sidediffs = diff(vertices);
        sidels = sqrt(sum(sidediffs.^2,2));
        tl = sum(sidels);
        sideps = round(sidels/tl*props.Density);
        boundary = [];
        for i=1:length(sidels)
            x = linspace(vertices(i,1),vertices(i+1,1),sideps(i));
            y = linspace(vertices(i,2),vertices(i+1,2),sideps(i));
            boundary = [boundary ; x' y'];
        end
        boundary = {boundary};
    case 'annulus'
        diameters = sort(pipeargs);
        d1 = diameters(1);
        d2 = diameters(2);
        th = linspace(0,2*pi,props.Density/2);
        th = th';
        boundary = {d1/2*[cos(th), sin(th)],d2/2*[cos(th), -sin(th)]};
    case 'wedge'
        assert(length(pipeargs)==2,'pipeCond:wedgeErr',...
            'Wedge arguments should be a list of radius, angle.');
        r = pipeargs(1);
        angle = pipeargs(2);
        assert(angle<2*pi,'pipeCond:wedgeErr',...
            'Wedge angle %f is not between 0 and 2pi');
        rpoints = round(props.Density/(2+angle));
        diampoints = round(angle*props.Density/(2+angle));
        x1 = linspace(0,r,rpoints);
        x2 = r*cos(linspace(0,angle,diampoints));
        x3 = linspace(x2(end),0,rpoints);
        y1 = zeros(1,rpoints);
        y2 = -r*sin(linspace(0,angle,diampoints));
        y3 = linspace(y2(end),0,rpoints);
        boundary = {[x1 x2 x3 ; y1 y2 y3]'};
    case 'boundary'
        boundary = pipeargs;
    case 'image'
        bwimage = pipeargs{1};
        spacing = pipeargs{2};
        boundary = getPerimeter(spacing,bwimage);
        tl = 0;
        for i=1:length(boundary)
            tl = tl + size(boundary{i},1);
        end
        ratio = props.Density/tl;
        if ratio < 2/3
            for i=1:length(boundary)
                a = boundary{i};
                l = size(a,1);
                np = l*props.Density/tl;
                skip = round(l/np);
                skip = skip + ~skip;
                aa = a(1:skip:l,:);
                boundary{i} = aa;
            end
        elseif ratio > 2
            numper = round(ratio);
            for i=1:length(boundary)
                a = boundary{i};
                aa = [];
                a = [a ; a(1,:)];
                for j=1:size(a,1)-1
                    x1 = a(j,1);
                    y1 = a(j,2);
                    x2 = a(j+1,1);
                    y2 = a(j+1,2);
                    xa = linspace(x1,x2,numper+1); xa = xa(1:end-1)';
                    ya = linspace(y1,y2,numper+1); ya = ya(1:end-1)';
                    aa = [aa ; [ xa ya ] ];
                end
                boundary{i} = aa;
            end
        end
    otherwise
        error('pipeCond:modeErr','pipeType ''%s'' not recognized.',pipeType)
end

%% Adjust length units
lengthUnits = struct('m',100,'in',2.54,'ft',30.48,'mm',0.1,'yd',91.44,...
    'km',1e5,'mi',1.609e5,'au',14.96e12,'ltyr',9.461e17,'hb',1.362e28,...
    'um',1e-4,'nm',1e-7,'cm',1);

if isfield(lengthUnits,props.Units)
    mult = lengthUnits.(props.Units);
else
    error('pipeCond:unitErr','Unit ''%s'' not recognized.',props.Units);
end


for i=1:length(boundary)
    boundary{i} = boundary{i}*mult;
end
L = L * mult;

%% Actual Area Inertia Call
I = getAreaInertia(boundary,props.ViewImage);

%% Get Gas and Velocity Info
if ~isnan(props.Gas)
    if ~isnan(props.Velocity)
        warning('pipeCond:gasWarn',['Gas specification ''%s'' will be' ...
            'overridden by Velocity specification.'],props.Gas);
    end
    
    gas2mass = struct('H2',2.016,'He',4.02,'N2',28.02,'O2',32,...
    'O3',48,'Ne',20.179,'Ar',39.948,'Kr',83.798,'Xe',131.293,'Rn',222,...
    'F2',38,'Cl2',70.91,'Br',79.9,'NH3',17.03,'C2H2',26,'air',29,...
    'C6H6',78.11,'CO2',44.01,'C4H10',58.1,'C2H6',30.07,'CH4',16.043,...
    'NO',30,'C3H8',44.09,'H2O',18.016,'Air',29,'C4H8',56.11,'CO',28.01,...
    'C2H4',28.03,'C6H12',84.16,'C6H14',86.17,'C7H16',100.2,...
    'CH3OH',32.04,'CH3Cl',50.49,'H2S',34.08,'NO2',46.01,'N2O',44.01,...
    'NO3',62.01,'C5H12',72.15,'C3H6',42.1,'S',32.06,'SO2',64.06,...
    'SO3',80.06,'SO',48.06,'C7H8',92.14,'C8H18',114.22,'CS2',76.13,...
    'SF6',146.05,'UF6',352,'XeF6',245.28,'OH',17.01,'OD',18.01,...
    'YO',104.91,'SrF',106.62,'CaF',59.08,'HfF',197.49,'ThO',251.04,...
    'Li',6.94,'Na',22.989,'K',39.1,'Rb',85.47,'Cs',132.9,'Sr',87.62,...
    'Dy',162.5,'Lu',174.97,'Er',167.26,'Ca',40.08);
    
    if isfield(gas2mass,props.Gas)
        m = gas2mass.(props.Gas);
    else
        error('pipeCond:gasErr','Gas ''%s'' not known.',props.Gas);
    end

    props.Mass = m;
end

if isnan(props.Velocity)
    v = sqrt(1.3801e-23 * props.Temp / ( 2*pi*props.Mass * 1.67e-27 ) );
else
    v = props.Velocity;
end

% 100 converts to cm/s
v = 100 * v;

%% Conductance at last

% divide by 1000 to convert from cm^3 to L
c = 1/2 * I/L * v / 1000; 




end