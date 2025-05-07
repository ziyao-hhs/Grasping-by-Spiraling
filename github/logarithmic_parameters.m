%% logarithmic spiral parameters
b = 0.16; 
interval = pi/6;
kp = 1/exp(b*interval);
theta01 = 3.5 * pi;
theta02 = 0.5 * pi; 
theta1 = theta01 : -interval : theta02 - 2*pi;
theta2 = theta01 : -interval : theta02;
k = 0.5 + 0.5 * exp (-2 * b * pi);

%% fist section radius R
ya = exp(b * theta01) * sin(theta01);
yb = k * exp(b * theta01) * sin(theta01);
deltay = abs(ya - yb);
xn = -exp(b * (theta01 - interval)) * cos(theta01 - interval);
yn = exp(b * (theta01 - interval)) * sin(theta01 - interval);
degree = atan(xn / (yn - ya));
R = deltay * sin(degree);


%% define contour logarithmic spiral
a = 44 / R; 
x1 = - a * exp (b * theta1) .* cos(theta1);
y1 = a * exp (b * theta1) .* sin(theta1);
x2 = - k * a * exp (b * theta2) .* cos(theta2);
y2 = k * a * exp (b * theta2) .* sin(theta2);

%% move spirals to origin
y_shift = y1(1); 
y1 = y1 - y_shift;
y2 = y2 - y_shift;

z1 = zeros(length(theta1));
z2 = zeros(length(theta2));

figure();
for i = 1:length(theta1)
    scatter3(x1(i), y1(i), z1(i), "green");
    hold on;
end

for i = 1:length(theta2)
    scatter3(x2(i), y2(i), z2(i), "red");
    hold on;
end

gamma = atan( y1(2) / x1(2));

%% define first section
p1 = [ x2(1), y2(1) - y2(1), z2(1)];
p2 = [ x1(13), y1(13) - y2(1), z1(13)];
p3 = [ x1(14), y1(14) - y2(1), z1(14)];
p4 = [ x2(2), y2(2) - y2(1), z2(2)];

p1_ = p1;
p2_ = p2 * rotate_z(gamma);
p3_ = p3 * rotate_z(gamma);
p4_ = p4 * rotate_z(gamma);
plot3([p1_(1), p2_(1),p3_(1),p4_(1),p1_(1)], [p1_(2), p2_(2),p3_(2),p4_(2),p1_(2)], [0,0,0,0,0]);

hold on;
scatter3(p1_(1),p1_(2), 0, 'green');hold on;
scatter3(p2_(1),p2_(2), 0, 'yellow');hold on;
scatter3(p3_(1),p3_(2), 0, 'red');hold on;
scatter3(p4_(1),p4_(2), 0, 'blue');hold on;

% the position of the threading hole
site1 = site(p1_, p2_, 5);
scatter3(site1(1),site1(2), 0, 'red');hold on;
site2 = site(p4_, p3_, 5);
scatter3(site2(1),site2(2), 0, 'blue');hold on;

site3 = site1 * rotate_x(deg2rad(120));
scatter3(site3(1),site3(2), site3(3), 'green');hold on;
site4 =site2 * rotate_x(deg2rad(120));
scatter3(site4(1),site4(2), site4(3),'yellow');hold on;

site5 = site1 * rotate_x(deg2rad(240));
scatter3(site5(1),site5(2), site5(3), 'black');hold on;
site6 =site2 * rotate_x(deg2rad(240));
scatter3(site6(1),site6(2), site6(3),'cyan');hold on;

xlabel('x');
ylabel('y');
zlabel('z');
h = gca;
plot3(h.XLim, [0 0], [0 0], 'r')
plot3([0, 0], h.YLim, [0 0], 'r');
plot3([0, 0], [0 0], h.ZLim, 'r');

% contour
plot( x1, y1, Color= "blue");
hold on;
plot( x2, y2, Color="red");

axis equal;

%% the total length L of the backbone
f = @(theta) sqrt((k * a * b * exp (b * theta)) .^ 2 + (k * a * exp (b * theta)) .^2);
L = integral(f, theta02, theta01);

%% the width D for the trunk tip
xDa = - a * exp (b * theta02) * cos(theta02);
yDa = a * exp(b * theta02) * sin(theta02);
xDb = - k * a * exp (b * theta02) * cos(theta02);
yDb = a * k * exp(b * theta02) * sin(theta02);
D = 2 * sqrt((xDa - xDb)^2 + (yDa - yDb)^2);

%% the width d of the smallest graspable range of the trunk body
xup = - a * exp(b * (theta02 - 2 * pi)) * cos(theta02 - 2 * pi);
yup = a * exp(b * (theta02 - 2 * pi)) * sin(theta02 - 2 * pi);
xdown = - a * exp(b * (theta02 - pi)) * cos(theta02 - pi);
ydown = a * exp(b * (theta02 - pi)) * sin(theta02 - pi);
d = sqrt((xup - xdown)^2 + (yup - ydown)^2);

%% the relationship between d/D
dDb = 1 / (exp(pi * b) - 1);


% rotation around x axis
function T_x = rotate_x(angle)
    T_x = [
    [1,0 0]
    [0, cos(angle), -sin(angle)]
    [0, sin(angle), cos(angle)]
    ];
end

% rotation around z axis
function T_z = rotate_z(angle)
    T_z = [
    [cos(angle),  -sin(angle), 0]
    [sin(angle), cos(angle), 0]
    [0, 0, 1]
    ];
end

function site_point = site(p1, p2, scale)
    site_point = [0,0,0];
    site_point(1) = p1(1) + (scale - 2) * ((p2(1) - p1(1)) / scale);
    site_point(2) = p1(2) + (scale - 2) * ((p2(2) - p1(2)) / scale);
end
