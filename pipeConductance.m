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
%    circle      radius
%
%    ellipse     [r1, r2]         either order
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
% 'Method'  Choose a specific method to use in calculating conductance. 
%           Can be 'Dushman', 'Knudsen', 'Numeric', or 'Clausing'. The latter
%           is only available for circular pipes, for which it is the
%           default method. Otherwise, 'Dushman' is default. Dushman
%           introduced the notion of considering short pipes as apertures
%           in series with a pipe. The result is exact as pipe shrinks to
%           zero or grows to infinity, and is a reasonable approximation in
%           between. You can specify 'Knudsen' to force the infinite pipe
%           approximation, but keep in mind that this always overestimates
%           conductance, even for rather long pipes. See Lafferty's
%           "Foundations of Vacuum Science and Technology", p 88.
%
%           If 'Numeric' is specified, molecules will be randomly
%           initialzied and propagated ti directly calculate the
%           transmission probability. This can be a slow process and poorly
%           converged for long pipes with low probabilities, in which case
%           'Dushman' is a much better choice.
%
% 'Number'  The number of molecules to initialize for the 'Numeric' Method,
%           see above. Default is 10,000.
%
% Outputs:
% c -       Conductance of the pipe in the rarified gas regime (i.e. dilute
%           molecular flow) measured in Liters / second.
%
% -------------------------------------------------------------------------
% EXAMPLE 1
% ---------
%   C = pipeConductance(1,'circle',1)
%   % Compare with formula from Pfieffer Vacuum:
%   % https://www.pfeiffer-vacuum.com/en/know-how/introduction-to-vacuum-technology/fundamentals/conductance/
%   % d = 1; l = 1; Cp = 12.1 * d^3 / l
%
% EXAMPLE 2
% ---------
%   C = pipeConductance('square',1,'ViewImage',true,'Density',60);
%   % 
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
    'Units','cm',...
    'Method','Dushman',...
    'Number',1e4);

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

area = nan;

assert(ischar(pipetype),'pipetype must be a string')
switch(pipetype)
    case 'circle'
        radius = pipeargs;
        th = linspace(0,2*pi,props.Density);
        th = th';
        boundary = {radius*[cos(th), -sin(th)]};
        area = pi*radius^2;
    case 'ellipse'
        r1 = pipeargs(1); r2 = pipeargs(2);
        th = linspace(0,2*pi,props.Density);
        th = th';
        boundary = {[r1*cos(th), -r2*sin(th)]};
        area = pi*r1*r2;
    case 'square'
        sidelength = pipeargs;
        a = round(props.Density/4);
        row = linspace(0,sidelength,a)';
        xp = [zeros(a,1) ; row ; ones(a,1)*sidelength ; flipud(row)];
        yp = [row ; ones(a,1)*sidelength ; flipud(row) ; zeros(a,1)];
        boundary = {[xp yp]};
        area = sidelength^2;
    case 'rectangle'
        pipeargs = sort(pipeargs);
        side1 = pipeargs(1);
        side2 = pipeargs(2);
        a = round(side1/(side1+side2)*props.Density/2);
        b = round(side2/(side1+side2)*props.Density/2);
        row = linspace(0,side1,a)';
        col = linspace(0,side2,b)';
        xp = [zeros(b,1) ; row ; ones(b,1)*side1 ; row(end:-1:1)];
        yp = [col ; ones(a,1)*side2 ; col(end:-1:1) ; zeros(a,1)];
        boundary = {[xp yp]};
        area = side1*side2;
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
        area = .5*side1*side2*sin(angle12);
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
        radii = sort(pipeargs);
        r1 = radii(1);
        r2 = radii(2);
        th = linspace(0,2*pi,props.Density/2);
        th = th';
        boundary = {r1*[cos(th), sin(th)],r2*[cos(th), -sin(th)]};
        area = pi*(r2^2-r1^2);
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
        area = 0.5*angle*r^2;
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

%% Prepare Area if not known already:
% Note Area is needed for the Dushman short length correction used by
% default.
if isnan(area)
    area = getArea(boundary);
end


%% Adjust units
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
area = area * mult * mult;



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
    v = sqrt(1.38065e-23 * props.Temp / ( 2*pi*props.Mass * 1.67262e-27 ) );
else
    v = props.Velocity;
end

% 100 converts to cm/s
v = 100 * v;



%% Conductance at last

switch(props.Method)
    case 'Knudsen'
        I = getAreaInertia(boundary,props.ViewImage);
        c = 1/2 * I/L * v / 1000; 
    case 'Dushman'
        I = getAreaInertia(boundary,props.ViewImage);
        c = 1/2 * I/L * v; 
        
        % Convert to unitless transition probability:
        alpha = c / area / v;
        
        % Dushman's aperture conductance
        alphad = 1/(1 + 1/alpha);
        c = alphad * area * v;
        
        % divide by 1000 to convert from cm^3 to L
        c = c / 1000 ;
    case 'Clausing'
        assert(strcmp(pipetype,'circle'),...
            'The Clausing method only applies to circular pipes.')
        f = @(r,x) x*sqrt(x^2+4*r^2)-x^2;
        g = @(r,x) (x^2-(2*x-1)*(x^2+4*r^2))/sqrt(x^2+4*r^2);
        aa = @(r,x) (f(r,1-x)-f(r,x))/(g(r,x)-g(r,1-x));
        r = radius; s7 = sqrt(7);
        a = aa(r/L,2*r*s7/(3*L+2*r*s7));
        q = sqrt(L^2+4*r^2);
        alpha = (1-2*a)/(3*r^2*L)*(4*r^3+(L^2-2*r^2)*q-L^3)...
            + a + (1-a)/(2*r^2)*(L^2-L*q+2*r^2);
        c = alpha * v * area / 1000;
    case 'Numeric'
        alpha = getTransmissionProb( boundary , L , props.Number);
        c = alpha * v * area / 1000;
        
end % end switch(props.Method)        


end % end pipeConductance function




%% Subfunctions Stored Internally for Convenient Packaging.

%% getPerimeter 
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

[B,~,N,~] = bwboundaries(crossArea,4);

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
        
    B{k} = a;

end

end

function A = getArea( Boundary )
    basex = Boundary{1}(1,1);
    basey = Boundary{1}(1,2);

    A = zeros(length(Boundary),1);
    
    for ii = 1:length(Boundary)
        points = Boundary{ii};
        assert(size(points,1)>1,'A Boundary needs at least two points')
        is2 = size(points,2);
        assert(is2==2,'Member #%d of Boundary has size Nx%d, not Nx2.',ii,is2)
        baserep = repmat([basex basey],size(points,1),1);
        pshift = points - baserep;
        chlengths = sqrt(sum(pshift.^2,2));
        chlengths2 = chlengths([2:end 1]);
        changles = atan2(pshift(:,2),pshift(:,1));
        dchangles = -diff(changles([1:end 1]));
        aaaa = .5.*chlengths.*chlengths2.*sin(dchangles);
        if ii==1
            A(ii) = sum(aaaa(2:end-1));
        else
            A(ii) = sum(aaaa);
        end
    end
    A = sum(A);
end

%% Boundary Sensibility Check
% Determines convexity, translates to first quadrant, removes doubled
% points, checks handedness and winding, checks array sizing.
function [boundary, convex] = checkBound( boundary )

    ll = length(boundary);
    assert(ll>0,'Boundary arraylist is empty.');

    xm = inf; 
    ym = inf;
    convex = false;

    ll = length(boundary);
    for i=1:ll

        % get the boundary in 'a', make sure it has the right size.
        a = boundary{i};
        istwo = size(a,2);
        assert(istwo==2,'Member #%d of Boundary has size Nx%d, not Nx2.',i,istwo)
        assert(size(a,1)>2,['Your Boundary needs at least three points' ...
            ' to enclose an area, and more will improve accuracy.']);


        % remove doubled points:
        aa = [diff(a); a(end,:)-a(1,:)];
        aa = abs(aa) > mean(abs(aa(:)))/100000;
        aa = aa(:,1) | aa(:,2);
        a = a(aa,:);


        % keep track of the minimum coordinates, so we can shift the boundary
        % to the first quadrant to remove any false zero coordinates.
        xm = min(xm,min(a(:,1)));
        ym = min(ym,min(a(:,2)));

        % Let's sum the exterior angles along the boundary to make sure it
        % doesn't overwrap, which would imply self intersection.
        aw = [a ; a(1,:)];
        sidediffs = diff(aw);
        sideangles = atan2(sidediffs(:,2),sidediffs(:,1));
        sideangles = [sideangles ; sideangles(1)];
        angles = mod(1e-6-diff(sideangles),2*pi);
        if ll==1
            convex = all(angles < pi);
        end
        angles = angles.*(angles < pi) - (2*pi - angles).*(angles > pi);
        winding = round(sum(angles)/(2*pi));
        if abs(winding) > 1
            error('areaInertiaErr:winding',['Sum of exterior angles'...
                ' along boundary is %dpi, indicating a self-overlapping'...
                ' and thus invalid boundary.'],2*winding)
        elseif winding < 0 && ll==1
            warning('areaInertiaWarn:winding',['Sum of exterior angles'...
                ' along boundary is %dpi, indicating that the boundary'...
                ' wraps counterclockwise. Reversing automatically, '...
                'but use caution with multi-boundaried areas.'],2*winding)
            a = flipud(a);
            convex = all(-angles < pi);
        end
        boundary{i} = a;
    end

    % On first loop we found the minimum coordinates, now we loop through and
    % translate all boundaries.
    offset = 1-[xm ym];
    for i=1:ll
        a = boundary{i};
        a = a + repmat(offset,size(a,1),1);
        boundary{i} = a;
    end
end %end checkbound subfunction

%% Prepare matrices representing chords of boundary
function chordData = prepChordMatrices( boundary, convex )
    % Perform a few preparatory tasks that involve looping through all the
    % boundaries.
    ll = length(boundary);
    lengths = zeros(ll,1);
    for i=1:ll

        % Populate the list of boundary lengths
        lengths(i) = size(boundary{i},1);
    end
    endps = cumsum(lengths);

    % Now I will populate eight matrices. The four matrices without 2 on the
    % end, together will represent all possible chords connecting two points in
    % the boundary. (A chord is specified by four points, hence the four
    % matrices). The matrices with 2's will also represent all possible chords,
    % but with a permutation of indices, which proves useful for easily
    % referencing a chord's closest neighbor later on. The particular
    % arrangement of these matrices of chords is important, and will be
    % discussed in further comments.
    xrep = zeros(endps(end),1);
    yrep = xrep;
    xrep2 = zeros(endps(end),1);
    yrep2 = xrep2;
    xtris = cell(1,ll);
    ytris = xtris;
    xtri2s = cell(1,ll);
    ytri2s = xtri2s;
    trackers = cell(1,ll);


    % for each connected sub-boundary in Boundary:
    for i=1:ll

        % get the boundary in 'a', make sure it has the right size.
        a = boundary{i};

        % get the x and y points along the boundary, make sure there are
        % enough. Translate the entire boundary into the first quadrant.
        x = a(:,1);
        y = a(:,2);
        l = length(x);

        % The 'rep' matrices will just be a column vector of boundary points
        % repeated across each row. Here we load the relevant slot with this
        % boundary.
        xrep(endps(i)-lengths(i)+1:endps(i)) = x;
        yrep(endps(i)-lengths(i)+1:endps(i)) = y;
        xrep2(endps(i)-lengths(i)+1:endps(i)) = x([end 1:end-1]);
        yrep2(endps(i)-lengths(i)+1:endps(i)) = y([end 1:end-1]);

        % Rearrange the perimeter points into matrices of points so that the i,jth
        % matrix element is the x (y) coordinate of the point j-1 to the right of
        % point i:
        % x1 x2 x3 x4 ...
        % x2 x3 x4 x1 ...
        % x3 x4 x1 x2 ...
        % ...
        xtri = reshape(repmat(x,l+1,1),l+1,l);
        ytri = reshape(repmat(y,l+1,1),l+1,l);
        xtri = xtri(1:end-1,:);
        ytri = ytri(1:end-1,:);

        % The tri2's are just like the tri's but with a unit permutation of the
        % x,y lists.
        xtri2 = reshape(repmat(x([end 1:end-1]),l+1,1),l+1,l);
        ytri2 = reshape(repmat(y([end 1:end-1]),l+1,1),l+1,l);
        xtri2 = xtri2(1:end-1,:);
        ytri2 = ytri2(1:end-1,:);
        
        % The tracker keeps track of the ordering of these matrices.
        tracker = endps(i)-lengths(i)+1:endps(i);
        tracker = tracker';
        tracker = reshape(repmat(tracker,l+1,1),l+1,l);
        tracker = tracker(1:end-1,:);

        % Load the tris into an arraylist. We'll build them into a larger 
        % matrix later.
        xtris{i} = xtri;
        ytris{i} = ytri;
        xtri2s{i} = xtri2;
        ytri2s{i} = ytri2;
        trackers{i} = tracker;
    end

    % Total length over all components:
    l = length(xrep);

    % Repeat rep matrices to their right:
    xrep = repmat(xrep,1,l);
    yrep = repmat(yrep,1,l);
    xrep2 = repmat(xrep2,1,l);
    yrep2 = repmat(yrep2,1,l);

    % Build larger matrices out of the xtris. This has the effect of making a
    % matrix like so:
    % a1 a2 a3 | b1 b2 b3 b4
    % a2 a3 a1 | b1 b2 b3 b4
    % a3 a1 a2 | b1 b2 b3 b4
    % --------------------
    % a1 a2 a3 | b1 b2 b3 b4
    % a1 a2 a3 | b2 b3 b4 b1
    % a1 a2 a3 | b3 b4 b1 b2
    % a1 a2 a3 | b4 b1 b2 b3
    %
    % Here ai's represent the ith x or y coordinates along the first boundary
    % in Boundary, bi's represent the ith x or y coordinates in the second
    % boundary, etc.
    xtri = blkdiag(xtris{:});
    ytri = blkdiag(ytris{:});
    xtri = xrep'.*(~xtri) + xtri;
    ytri = yrep'.*(~ytri) + ytri;
    xtri2 = blkdiag(xtri2s{:});
    ytri2 = blkdiag(ytri2s{:});
    xtri2 = xrep2'.*(~xtri2) + xtri2;
    ytri2 = yrep2'.*(~ytri2) + ytri2;
    tracker = blkdiag(trackers{:});
    tracker = repmat(1:l,l,1).*(~tracker)+tracker;
    
    % Now we want to rearrange the matrix like so:
    % a1 a2 a3 | b1 b2 b3 b4 | c1 c2
    % a2 a3 a1 | b1 b2 b3 b4 | c1 c2
    % a3 a1 a2 | b1 b2 b3 b4 | c1 c2
    % ----------------------------
    % b1 b2 b3 b4 | a1 a2 a3 | c1 c2
    % b2 b3 b4 b1 | a1 a2 a3 | c1 c2
    % b3 b4 b1 b2 | a1 a2 a3 | c1 c2
    % b4 b1 b2 b3 | a1 a2 a3 | c1 c2
    % ------------------------------
    % c1 c2 | a1 a2 a3 | b1 b2 b3 b4
    % c2 c1 | a1 a2 a3 | b1 b2 b3 b4
    %
    % That is to say, we want to move the nth block from the diagonal, where n
    % is the nth boundary in Boundary, adjacent to the left edge of the matrix,
    % and slide everything else over in its wake. The following achieves this:
    for i=1:ll
        e = endps(i);
        b = e - lengths(i)+1;
        xtri(b:e,1:e) = xtri(b:e,[b:e 1:b-1]);
        ytri(b:e,1:e) = ytri(b:e,[b:e 1:b-1]);
        xtri2(b:e,1:e) = xtri2(b:e,[b:e 1:b-1]);
        ytri2(b:e,1:e) = ytri2(b:e,[b:e 1:b-1]);
    end
    
    % Why did we do all that? The four matrices xrep, yrep, xtri, and ytri now
    % represent every possible chord connecting a pair of points in
    % Boundary. For a given pair (i,j), the chord connecting the points
    % (xrep(i,j),yrep(i,j)) and (xtri(i,j),ytri(i,j)) connects the ith point in
    % the Boundary, where i steps through each boundary in Boundary in a
    % clockwise ordered manner, to the point j away, where j steps first
    % through all the points in the same sub-boundary that i is in, beginning
    % the next one clockwise from i, and then through all points in the other
    % sub-boundaries, visiting the sub-boundaries at random but stepping
    % clockwise within each sub-boundary.
    %
    % The 2's are the same, but with a coordinate permutation.

    % We can now get the length of these chords easily:
    chlen = sqrt((xtri - xrep).^2 + (ytri - yrep).^2);

    % We can also get the angle they make with the x-axis:
    chang = atan2(ytri-yrep,xtri-xrep);
    chang2 = atan2(ytri2-yrep2,xtri2-xrep2);
    chang3 = atan2(ytri2-yrep,xtri2-xrep);

    % Since the first column of chords connect each point to itself, they have 
    % length zero. The second column represent neighboring chords, which are
    % not really chords but the segments making up the perimeter. The angles in
    % the second column of chang therefore represent the angle of the
    % perimeter with the x-axis. This will be used to find the cosine of the
    % angle of each chord with the normal to the perimeter in the integrand we
    % seek.
    changt = repmat(chang(:,2),1,l) - chang;
    changt2 = repmat(chang2(:,2),1,l) - chang;
    changt = mod(changt,2*pi);
    changt2 = mod(changt2,2*pi);

    % Since the integrand proceeds over angles emanating from each point along
    % the perimeter, we'll want to weight each chord by the angular
    % displacement to the next chord (i.e. a right Reimann sum). We could get
    % this by computing a diff over adjacent columns of the chord angles with
    % the x-axis. Problems would occur along the boundaries of the master chord
    % matrix where chords switch between different sub-boundaries. To address
    % this, we don't use diff, but instead compare relative to the permuted
    % matrix of chords. This ensures that each chord is compared to a chord
    % emanating from the same point and connected to the adjacent point on the
    % same sub-boundary.
    changd = chang - chang3;

    % Now changd should be small but of either sign, and there could be
    % wraparound errors if chords and their neighbors are on different sides of
    % horizontal. We want to make sure we always think of it as small and
    % positive:
    changd = mod(changd,2*pi);
    changd = changd.*(changd<pi) + (2*pi-changd).*(changd>pi);

    if ~convex

        % Conveniently, we can reject one class of chords that are outside
        % the area of interest purely based on whether its angle is between the
        % two boundary edges that it shares a vertex with:
        chout1 = changt <= pi & changt2 <= pi; %  chord included if true.

        % Lets also note which chords intersect an internal exclusion. We do this
        % by checking for intersections with the perimeter. This is the longest
        % aspect of the calculation by far, and so a waitbar is implemented. I
        % thought of implementing this with MATLAB's 'inpolygon' function, but it
        % proved a good bit slower than the following. I suspect it may be quicker
        % to choose points along each chord and see if any 'hit' a boolean matrix
        % encoding the invalid region as an area. This should be O(n^2) instead of
        % O(n^3), but I don't know how the overhead compares.
        h = waitbar(0, 'Rejecting Invalid Chords ...');
        chout2 = true(l,l);
        % for each perimeter segment:
        for i=1:l

            % The line defined by (xs1,ys1),(xs2,ys2) is a segment in the boundary.
            xs1 = xrep(i,2);
            ys1 = yrep(i,2);
            xs2 = xtri(i,2);
            ys2 = ytri(i,2);

            % The following checks the sign of the cross product of the vector from
            % xs1,ys1 to both endpoints of the chord. If they are opposite, the
            % endpoints are on opposite sides of the segment, so the chord crosses 
            % this perimeter segment.
            chord_x_segment = sign((xs2-xs1)*(yrep-ys1)-(xrep-xs1)*(ys2-ys1))==-sign((xs2-xs1)*(ytri-ys1)-(xtri-xs1)*(ys2-ys1));

            % The following checks whether the two endpoints of this perimeter
            % segment also lie on opposite sides of the chord. The '~=' means we
            % also include the case where one of the endpoints of this perimeter
            % segment is on top of the chord, since when this happens the sign
            % argument and output are zero. We didn't include this case above
            % because we don't want to exclude chords that share an endpoint with
            % perimeter segments (since all of them do that!).
            segment_x_chord = sign((xtri-xrep).*(ys1-yrep)-(ytri-yrep).*(xs1-xrep))~=sign((xtri-xrep).*(ys2-yrep)-(ytri-yrep).*(xs2-xrep));

            % if each crosses the other, we have an intersection, so we reject the
            % chord by flipping the bit corresponding to it in the chout2.
            intx = chord_x_segment & segment_x_chord;
            chout2 = chout2 & ~intx;

            % Update the waitbar
            waitbar( i / l, h);
        end
        close(h)

    else
        chout1 = 1;
        chout2 = 1;
    end
    
    chordData.include = chout1.*chout2;
    chordData.length = chlen;
    chordData.tangent = changt;
    chordData.dangle = changd;
    chordData.x1 = xrep(:,2);
    chordData.x2 = xtri(:,2);
    chordData.y1 = yrep(:,2);
    chordData.y2 = ytri(:,2);
    chordData.track = tracker;

end


%% getAreaInertia
function I = getAreaInertia( Boundary, verbose )
%% This function computes the "area inertia" of a simply connected but
% otherwise arbitrary pipe cross section. I invented the term area inertia,
% but you can find more information in W. Steckelmacher's 1966 review:
% "A review of the molecular flow conductance for systems of tubes and 
% components and the measurement of pumping speed"
% See page 5, equation 14. (which has an error by the way. chord length
% should be squared)
%
% Boundary is an arraylist of Nx2 arrays, each of which is an ordered list
% of points along a single connected component of the boundary, with the
% area on the right (i.e. oriented clockwise). If only a single component,
% counterclockwise orientation is detected and corrected automatically.
%
% verbose produces a figure if true. This is best used with a total number
% of boundary vertices less than 100, or the figure will take too long.

[Boundary, convex] = checkBound( Boundary );
cc = prepChordMatrices( Boundary, convex );


% spacing holds the length of each perimeter segment. This is "ds" in
% an integral.
space = cc.length(:,2);
space = repmat(space,1,length(space));


% integral over perimeter and chord angle of half chord length squared
% cosine of the angle to the normal (sine of angle to perimeter tangent)
% with non-internal chords excluded:
integrand = .5*cc.length.^2.*sin(cc.tangent).*cc.dangle.*space.*cc.include;

% We reject the first two columns of chords, since they are of length zero
% and on the boundary, respectively.
I = sum(sum(integrand(:,1:end)));


% The following is a slow visual check that the right chords are included.
if verbose
    figure;
    hold on;
    for i=1:l
        for j=1:l
            if chout1(i,j) && chout2(i,j)
                plot([xrep(i,j) xtri(i,j)],[yrep(i,j) ytri(i,j)])
                pause(.01)
            end
        end
    end
end

% This is a slow check to make sure that all chord neighbors are chosen
% properly. It focuses on chord neighbors more than pi/10 degrees away.
if false
    figure;
    hold on;
    for i=1:ll
        plot(Boundary{i}(:,1),Boundary{i}(:,2))
    end
    for i=1:l
        for j=1:l
            if changd(i,j) > pi/10
                line([xrep(i,j) xtri(i,j)],[yrep(i,j) ytri(i,j)],'Color','b')
                hold on
                fill([xrep(i,j) xtri(i,j) xtri2(i,j)],[yrep(i,j) ytri(i,j) ytri2(i,j)],'r')
                hold off
                input('');
            end
        end
    end
end

end % end getAreaInertia sub function



function p = getTransmissionProb( Boundary, L, N )

    [Boundary, convex] = checkBound( Boundary );
    cc = prepChordMatrices( Boundary, convex );
    
    
    maxX = max(cc.x1);
    maxY = max(cc.y1);
    minX = min(cc.x1);
    minY = min(cc.y1);
    l = length(cc.x1);
    
    % Choose initial coordinates within the bounded area
    initx = [];
    inity = [];
    
    % How many to try at once:
    M = 1000;
    
    % Loop until a full N are obtained.
    while length(initx) < N
        % guess points anywhere within the smallest square that encloses
        % the boundary.
        guessx = rand(M,1)*(maxX-minX)+minX;
        guessy = rand(M,1)*(maxY-minY)+minY;
        
        % make M x B matrices, where B is the number of points on the
        % perimiter of the boundary.
        guessxM = repmat(guessx,1,l);
        guessyM = repmat(guessy,1,l);
        x1M = repmat(cc.x1',M,1);
        x2M = repmat(cc.x2',M,1);
        y1M = repmat(cc.y1',M,1);
        y2M = repmat(cc.y2',M,1);
        
        % for each point, see whether the line due north from that point
        % crosses the boundary perimeter an odd number of times.
        doesSpan = (x1M <= guessxM & guessxM < x2M) | ...
                   (x2M <= guessxM & guessxM < x1M) ;
        yeff = (guessxM - x1M)./(x2M - x1M).*(y2M-y1M)+y1M;
        isAbove = yeff > guessyM;
        numcross = sum(isAbove & doesSpan,2);
        keep = ~~mod(numcross,2);
        initx = [initx ; guessx(keep)];
        inity = [inity ; guessy(keep)];
    end
    initx = initx(1:N);
    inity = inity(1:N);
    
%     figure;hold on
%     plot(cc.x1,cc.y1,'b-')
%     plot(initx,inity,'r*','LineStyle','none')
    
    % Now we choose direction. Can be anything in the 2pi of solid angle
    % that points into the pipe, but steeper angles are less likely by
    % cos(th) because the cross section looks smaller for molecules coming
    % at those angles.
    th = acos(2*rand(N,1)-1)/2; % polar angle
    ph = rand(N,1)*2*pi;  % azimuthal angle
    
    % Find boundary contact in this direction:
    x1M = repmat(cc.x1',N,1);
    x2M = repmat(cc.x2',N,1);
    y1M = repmat(cc.y1',N,1);
    y2M = repmat(cc.y2',N,1);
    phM = repmat(ph,1,l);
    xxM = repmat(initx,1,l);
    yyM = repmat(inity,1,l);
    xpM = xxM + cos(phM);
    ypM = yyM + sin(phM);
    
    % Now we determine contact between the i'th molecule and the j'th
    % boundary segment and store it in the element (i,j) of an NxM boolean
    % matrix. The molecules position is (xxM,yyM), and another point in its
    % direction of motion is (xpM,ypM). (x1M,y1M) and (x2M,y2M) bound the
    % perimeter segment in question.
    xv = xpM - xxM;     yv = ypM - yyM;
    x1v = x1M - xxM;    y1v = y1M - yyM;
    x2v = x2M - xxM;    y2v = y2M - yyM;
    side1 = sign(xv.*y1v-yv.*x1v);
    side2 = sign(xv.*y2v-yv.*x2v);
    contact = side1 ~= side2;
    
    % Now we have to eliminate backwards contacts:
    distance = sqrt( x1v.^2 + y1v.^2);
    dist1 = distance+sqrt(x2v.^2+y2v.^2);
    x1vp = x1v-cos(phM)*1e-10;
    x2vp = x2v-cos(phM)*1e-10;
    y1vp = y1v-sin(phM)*1e-10;
    y2vp = y2v-sin(phM)*1e-10;
    dist2 = sqrt(x1vp.^2+y1vp.^2)+sqrt(x2vp.^2+y2vp.^2);
    backwards = dist1 < dist2;
        
    % Now we need the distance to said contact.
    distance(~contact) = inf;
    distance(backwards) = inf;
    [distlist, perpos] = min(distance,[],2);
    pipepos = cos(th).*distlist;
    
%     % Let's check it in a figure
%     figure; hold on
%     %plot(cc.x1,cc.y1,'k-')
%     for i=1:N
%         line([initx(i),initx(i)+cos(ph(i))/10],[inity(i),inity(i)+.1*sin(ph(i))])
%         plot(initx(i)+cos(ph(i))/10,inity(i)+sin(ph(i))/10,'b+')
%         plot(cc.x1(perpos(i)),cc.y1(perpos(i)),'r*')
%         pause(.01)
%     end
    
    % Now we do the loop
    valid_angles = cc.tangent;
    valid_angles(~cc.include) = nan;
    valid_angles(:,1) = nan;
    
    throughput = sum(pipepos >= L);
    perpos = perpos(pipepos < L);
    pipepos = pipepos(pipepos < L);
    M = length(pipepos);
    
    while ~isempty(pipepos)
        % choose next direction. This time must follow the lambert cosine
        % emission law. Since we're in 3D, there is still the sin(th) in
        % the jacobean, so really we want a probability distribution 
        % p(th) ~ sin(th)cos(th) ~ sin(2*th) on 0 to pi/2. 
        % The following achieves this, as you can confirm by plotting:
        % hist(acos(rand(1e6,1)*2-1)/2,100)
        % To derive it, integrate and invert the CDF of the desired density.
        th = acos(rand(M,1)*2-1)/2;
        ph = rand(M,1)*2*pi;
        
        angle_plain = acot(cos(th)./(cos(ph).*sin(th))) + pi/2;
        angle_pipes = asin(sin(ph).*sin(th));
        
        anglemat = valid_angles(perpos,:);
        
        [~,loc] = min((anglemat - repmat(angle_plain,1,l)).^2,[],2);

        proplength = cc.length(sub2ind([l l],perpos,loc));
        proplength = proplength.*tan(angle_pipes);
        pipepos = pipepos + proplength;
        
        perpos = cc.track(sub2ind([l l],perpos,loc));
        
        throughput = throughput + sum(pipepos >= L);
        perpos = perpos(pipepos < L & pipepos > 0);
        pipepos = pipepos(pipepos < L & pipepos > 0);
        M = length(pipepos);
    end
    

    p = throughput / N;

    
    % loop until enough valid points obtained:
    %  choose initial points somewhere in pipe cross sectional area
    % end loop
    % choose initial direction
    % find first boundary contact
    % loop:
    %  note if return or exit has occurred
    %  choose next direction
    % end loop

end