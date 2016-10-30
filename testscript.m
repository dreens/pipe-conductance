%% Here I compare with the circular formula:
ss = pi/100;
th = ss:ss:2*pi;
B = {1.4287*[cos(th') -sin(th')]}; 
I = getAreaInertia(B)
%R = getAreaInertia({10*B{1}})/I
getConductance(4.65e-26,293,10,I)
12.1*(2*1.4287)^3/10

%% And here I compare with the annulus formula:
ss = pi/100;
th = (ss:ss:2*pi)';
B = {[cos(th) sin(th)],[2*cos(th) -2*sin(th)]};
I = getAreaInertia(B)
getConductance(4.65e-26,293,10,I)
8*pi*1.154*sqrt(1.38e-23*293/(2*pi*4.65e-26))*.01
% seems to be off by a factor of 3*pi/8 for some reason. Close enough?

%% Do repeated points matter?
% No they don't.
ss = pi/100;
th = ss:ss:2*pi;
B = {[cos(th') -sin(th')]};
I = getAreaInertia(B)
th = [th(1:30) th(30)*ones(1,20) th(31:end)];
B = {[cos(th') -sin(th')]};
I = getAreaInertia(B)

%% Do disjoint areas give the conductance sum?
ss = pi/100;
th = ss:ss:2*pi;
B = {[cos(th') -sin(th')],[1.1+cos(th'), -2.5-sin(th')]};
I = getAreaInertia(B)
B = {[cos(th') -sin(th')]};
I = getAreaInertia(B)


%% Let's do the final calculation now.
% (.313,.313)->(.313,.847)
s = .28*2.54;
b = .847*2.54;
r = .25*2.54;
th = linspace(0,pi,50);
x1 = linspace(s,b,50);
y1 = ones(1,50)*s;
x2 = b+s*cos(pi/2-th);
y2 = s*sin(pi/2-th);
x3 = fliplr(x1);
y3 = -y1;
x4 = [y1 y2 y3];
y4 = -[x1 x2 x3];
x5 = -[x1 x2 x3 x4];
y5 = -[y1 y2 y3 y4];
x6 = r*cos([th, th+pi]);
y6 = r*sin([th, th+pi])+b;
x7 = x6+b;
y7 = y6-b;
x8 = -x6;
y8 = -y6;
x9 = -x7;
y9 = -y7;
x = [x1 x2 x3 x4 x5 x6 x7 x8 x9];
y = [y1 y2 y3 y4 y5 y6 y7 y8 y9];
figure;plot(x,y)
B = {[x1 x2 x3 x4 x5 ; y1 y2 y3 y4 y5]',[x6;y6]',[x7;y7]',[x8;y8]',[x9;y9]'};
I = getAreaInertia(B);
getConductance(4.65e-26,293,20,I)



