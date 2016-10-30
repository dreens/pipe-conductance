%% Check clockwise orientation
% It does seem to be clockwise. I posted on the file exchange to double
% check.
a = rand(401) > 0.992;
a = imdilate(a,strel('disk',5));
imtool(a)
b = bwboundaries(a);
figure;hold on;
for j=1:length(b)
    bb = b{j};
    for i = 1:length(bb)-1
        plot(bb(i:i+1,1),bb(i:i+1,2))
        pause(.02)
    end
end

%% See about smoothing a boundary
% Can be done, maybe later though.
a = fspecial('disk',100) > 0;
b = bwboundaries(a);
b = b{1};
xx = b(:,1);
yy = b(:,2);
cftool

%% Let's test getPerimeter
% Let's open the profile of interest. I obtained this by making a part in
% solidworks, viewing it straight on, coloring it low-gloss black plastic,
% setting an all-white scene background, choosing the highest resolution in
% System Options, Document Properties, Image Quality, and then saving as
% a PNG with black and white resolution and 1200dpi set in the PNG options
% button of the save dialog box.
aa = imread('R:\Data Analysis\Pipe Conductance\Matlab_Profile.PNG');
aa = aa==0;
% Now to check the spacing I count pixels in a half inch: 6429-4541. I did
% this by zooming in on a feature known to be 1/2" (1.27cm) with imtool.
%imtool(aa)
%cm_per_pixel_y = 1.27 / ( 6679 - 4448 );
cm_per_pixel_x = 1.27 / ( 5763 - 3874 );

b = getPerimeter(cm_per_pixel_x, cm_per_pixel_x, aa);
%% Plot the edges
% The handedness of each boundary is indicated by increasing the color
% along the boundary gradually starting from black.
figure;
colors = get(gca,'ColorOrder');
hold on
for j=1:length(b)
    for kj = 1:length(b{j})-1
        plot(b{j}(kj:kj+1,1),b{j}(kj:kj+1,2),'Color',colors(j,:)*kj/length(b{j}))
    end
end
axis equal

%% Some math by hand
% without pins, I find C = 8.1 over 20cm. With a pin pair, I find C = 2.6
% over the same distance. In practice, we can say the pins are ~2.5mm
% "thick" since you could replace a 3.175mm diameter circle with a 2.5 by
% 3.175mm circle and have the same area. The stages are about 5.5mm apart.
% This means that 2.5/5.5 = 45% of the tube has the more restrictive
% conductance and 55% has the lesser. Taking the weighted harmonic mean:
% C = 100/( 45/2.6 + 55/8 ) = 4.1 L / s.
%
% By the way, without pins, conductance would be 28.2 for a circular tube
% with the same cross sectional area.

%% Pass it on and get Conductance.
I = getAreaInertia(b)
C = getConductance(4.65e-26,293,20,I)