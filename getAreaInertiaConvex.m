function I = getAreaInertiaConvex( xp, yp )
% This function computes the "area inertia" of a convex but otherwise
% arbitrary pipe cross section. I invented the term "area inertia", but you
% can find more information in W. Steckelmacher's 1966 review titled:
% "A review of the molecular flow conductance for systems of tubes and 
% components and the measurement of pumping speed"
% See page 5, equation 14.

% Make column vectors
xp = xp(:);
yp = yp(:);

% Length is always useful
l = length(xp);

% A matrix of chord lengths between each point pair around the perimeter:
xpb = repmat(xp,1,l); ypb = repmat(yp,1,l);
chordlength = sqrt((xpb-xpb').^2+(ypb-ypb').^2);

% and the special case of the perimeter:
spacing = sqrt(diff([xp; xp(1)]).^2 + diff([yp; yp(1)]).^2);
spacing = repmat(spacing,1,l-2);

% lets rearrange chordlength so that i,j has the distance from the ith
% point to the point j to the right of it.
for j=1:l
    chordlength(j,:) = chordlength(j,[j:end 1:j-1]);
end

% Now we'll need the angles between chords emanating from each perimeter
% point. To get this, we'll get the sidelengths of triangles enclosing the
% desired angles and invert the law of cosines.
tri1 = chordlength(:,2:(end-1));
tri2 = chordlength(:,3:end);
tri3 = tri2;

mod1 = @(x,y)mod(x-1,y)+1;
for i=1:l
    for j=1:l-2
        tri3(i,j) = chordlength(mod1(i+j,l),2);
    end
end

angle = acos((-tri3.^2+tri2.^2+tri1.^2)./(2*tri1.*tri2));

% Here if I could somehow zero angles corresponding to concave chord
% sections this could generalize to concave areas...

cangle = cumsum(angle,2);

integral = .5*chordlength(:,1:end-2).*sin(cangle).*angle.*spacing;
    
I = sum(sum(integral));
