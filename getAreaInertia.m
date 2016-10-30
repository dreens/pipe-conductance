function I = getAreaInertia( Boundary )
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
% area on the right (i.e. oriented clockwise). 

% Convert boundary into a master list of x and y coordinates:
ll = length(Boundary);
assert(ll>0,'Boundary arraylist is empty.');

lengths = zeros(ll,1);
for i=1:ll
    lengths(i) = size(Boundary{i},1);
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

% for each connected sub-boundary in Boundary:
for i=1:ll
    
    % get the boundary in 'a', make sure it has the right size.
    a = Boundary{i};
    istwo = size(a,2);
    assert(istwo==2,'Member #%d of Boundary has size Nx%d, not Nx2.',i,istwo)
    
    % get the x and y points along the boundary, make sure there are
    % enough.
    x = a(:,1);
    y = a(:,2);
    l = length(x);
    assert(l>2,['Your Boundary needs at least three points' ...
        ' to enclose an area, and more will give more accurate results.']);
    
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
    
    % Load the tris into an arraylist. We'll build them into a larger 
    % matrix later.
    xtris{i} = xtri;
    ytris{i} = ytri;
    xtri2s{i} = xtri2;
    ytri2s{i} = ytri2;
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

% Conveniently, we can reject one class of chords that are outside
% the area of interest purely based on whether its angle is between the
% two boundary edges that it shares a vertex with:
chout1 = changt < pi & changt2 < pi; %  chord included if true.

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

% spacing holds the length of each perimeter segment. This is "ds" in the
% integral.
spacing = chlen(:,2);
spacing = repmat(spacing,1,l);

% integral over perimeter and chord angle of half chord length squared
% cosine of the angle to the normal (sine of angle to perimeter tangent)
% with non-internal chords excluded:
integrand = .5*chlen.^2.*sin(changt).*changd.*spacing.*chout1.*chout2;

% We reject the first two columns of chords, since they are of length zero
% and on the boundary, respectively.
I = sum(sum(integrand(:,3:end-1)));


% The following is a slow visual check that the right chords are included.
if false
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
