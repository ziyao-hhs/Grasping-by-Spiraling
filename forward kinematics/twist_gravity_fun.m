function [third_line_change, new_point] = twist_gravity_fun(D, D2, gravity_direction)
b = 0.16; 
interval = pi/6; 
kp = 1/exp(b*interval);
theta01 = 3.5 * pi;
theta02 = 0.5 * pi +  interval;
theta1 = theta01 : -interval : theta02 - 2*pi;
theta2 = theta01 : -interval : theta02;
k = 0.5 + 0.5 * exp (-2 * b * pi);
section_number = 18;
alpha_value = 1; % visulization transparency level

%% the radius R for the first section 
ya = exp(b * theta01) * sin(theta01);
yb = k * exp(b * theta01) * sin(theta01);
deltay = abs(ya - yb);
xn = -exp(b * (theta01 - interval)) * cos(theta01 - interval);
yn = exp(b * (theta01 - interval)) * sin(theta01 - interval);
degree = atan(xn / (yn - ya));
R = deltay * sin(degree);


%% define contour logarithmic spiral
a = 44 / R; % make sure R = 44 mm
x1 = - a * exp (b * theta1) .* cos(theta1);
y1 = a * exp (b * theta1) .* sin(theta1);
z1 = zeros(1, length(theta1));

x2 = - k * a * exp (b * theta2) .* cos(theta2);
y2 = k * a * exp (b * theta2) .* sin(theta2);
z2 = zeros(1, length(theta2));

% length of the sections
L = zeros(1,section_number);
l1 = sqrt(x2(2)^2 + (y2(2) - y2(1))^2 + z2(2)^2);
sec1_length = sqrt(x2(2)^2 + (y2(2) - y2(1))^2 + z2(2)^2);
sitex = zeros(1,section_number);
sitex(1) = sec1_length;
L(1) = sec1_length;


% right side length of the first section
side1 = sqrt((x2(2) - x1(14))^2 + (y2(2) - y1(14))^2 + (z2(2) - z1(14))^2);

%% positions for the threading holes
sec1_length = 52.1271;
sec1_height = 44;
sec1_cable_y = sec1_height;
upper_cable_x = 19.4797;
lower_cable_x = 47.4093;
ori_coordinate = [[1, 0, 0, 0] 
                  [0, 1, 0, 0]
                  [0, 0, 1, 0] 
                  [0, 0, 0, 1]];

site1 = [upper_cable_x, sec1_cable_y, 0];
site2 = [upper_cable_x, -sec1_cable_y/2, sqrt(3) * sec1_cable_y / 2];
site3 = [upper_cable_x, -sec1_cable_y/2, - sqrt(3) * sec1_cable_y / 2];

X = [lower_cable_x, sec1_cable_y, 0];
Y = [lower_cable_x, -sec1_cable_y/2, sqrt(3) * sec1_cable_y / 2];
Z = [lower_cable_x, -sec1_cable_y/2, - sqrt(3) * sec1_cable_y / 2];

T1 = ori_coordinate * trans(sec1_length, 0, 0);
T1_point = [T1(1,4), T1(2,4), T1(3,4)];
site4 = 6/10 * X + 4/10 * T1_point;
site5 = 6/10 * Y + 4/10 * T1_point; 
site6 = 6/10 * Z + 4/10 * T1_point;

site1 = 6/10 * site1;
site2 = 6/10 * site2;
site3 = 6/10 * site3;

for i = 2:section_number
    L(i) = l1 * kp ^ i; 
    sitex(i) = sitex(i - 1) + l1 * kp ^ i;
end

% beam theory parameters
%% bending moment I (stiffness)
radius = zeros(1 , (section_number - 1));
r = 3.825 * 1e-3;
for i = 1 : (section_number - 1)
    radius(i) = r * (kp ^ (i - 1));
end

stiffness = zeros(1 , (section_number - 1));
for i = 1 : (section_number - 1)
    stiffness(i) = (pi / 4) * ( radius(i) ^ 4 );
end

%% Youngs Modulus E
E = 16.9e6;

%% moment M (Torque)
Torque = zeros(1 , (section_number - 1));
torque = 0.0127;
for i = 1 : (section_number - 1)
    Torque(section_number - i) = torque * (kp ^ (-1/2) ^ (i - 1));
end

%% joint length L 
L_joint = zeros(1 , (section_number - 1));
L_init = 2.10 * 1e-3;
for i = 1 : (section_number - 1)
    L_joint(i) = L_init * (kp ^ (i - 1));
end

%% changing angles: theta_bend, and end_length
theta_bend = zeros(1 , (section_number - 1));
for i = 1 : (section_number - 1)
    theta_bend(i) = Torque(i) * L_joint(i) / (E * stiffness(i));
end

r_hole = zeros(1 , (section_number - 1));
r_hole_init = 28.87 * 1e-3;
for i = 1 : (section_number - 1)
    r_hole(i) = r_hole_init * kp ^ i;
end

end_length = zeros(1 , (section_number - 1));
for i = 1 : (section_number - 1)
    end_length(i) = sqrt(2 * r_hole(i) ^ 2 * (1 - cos(pi / 6 + theta_bend(i))));
end

%% original cable length between sections (gravity = 90°)
delta_length = zeros(1 , (section_number - 1));
origin_length = zeros(1 , (section_number - 1));
origin_length_init = 13.74 * 1e-3;
for i = 1 : (section_number - 1)
    origin_length(i) = origin_length_init * kp ^ (i - 1);
end

for i = 1 : (section_number - 1)
    delta_length(i) = end_length(i) - origin_length(i);
end

%% portion
Portion = zeros(1 , (section_number - 1));
for i = 1 : (section_number - 1)
    Portion(i) = delta_length(i) / sum(delta_length);
end

hole_point = 0.6;
t = zeros(1, (section_number - 1));
a = hole_point * side1;
b = hole_point * side1;

t(section_number - 1) = sqrt(a ^ 2 + b ^ 2 - 2 * a * b * cos(pi/6));
t_sum = sum(t);

for i = (section_number - 2) : -1 : 1
    t(i) = t(i + 1) * kp;
    t_sum = t_sum + t(i);
end

%% intial cable length under gravity
gravity_angle = solve_angle(deg2rad(gravity_direction));
gravity_tendon_length = zeros(1,section_number - 1);
for i = 1 : section_number - 1
    a = r_hole(i);
    b = r_hole(i);
    gravity_tendon_length(i) = angle2tendon(a, b, pi/6 + deg2rad(gravity_angle(i)));
end

%% first cable rotation angle
new_tendon_length = gravity_tendon_length;
distribute_length = zeros(1, (section_number - 1));
extrem_flag = zeros(1 , (section_number - 1));
for i = (section_number - 1) : -1 : 1 
    distribute_length(i) = distribute_length(i) + D * 1e-3 * Portion(i);
    if distribute_length(i) < gravity_tendon_length(i)
        new_tendon_length(i) = new_tendon_length(i) - distribute_length(i);
    else 
        new_tendon_length(i) = 0;
        addtion = distribute_length(i) - gravity_tendon_length(i);
        sum_rest = sum(delta_length(1 : i - 1));
        extrem_flag(i) = 1;
        for j = 1 : i - 1
            distribute_length(j) = distribute_length(j) + addtion * delta_length(j) / sum_rest;
        end
    end
end

final_angle = zeros(1, (section_number - 1));
for k = 1 : (section_number - 1)
    a = r_hole(k);
    b = r_hole(k);
    final_angle(k) = tendon2angle(a, b, new_tendon_length(k));
end

rotate_angle = zeros(1, (section_number - 1));
for k = 1 : (section_number - 1)
    rotate_angle(k) = pi/6 - final_angle(k);
end

new_coordinate = zeros(section_number, 4, 4);
new_point = zeros(section_number, 4);

last_point1 = zeros(section_number - 1, 4);
last_point2 = zeros(section_number - 1, 4);
last_point3 = zeros(section_number - 1, 4);

active_point = zeros(section_number - 1, 4);
passive_point1 = zeros(section_number - 1, 4);
passive_point2 = zeros(section_number - 1, 4);

passive_tendon_length = zeros(section_number - 1,1);
active_tendon_length = zeros(section_number - 1,1);

%% actuating the first cable
for k = 1 : (section_number)
    if k == 1 
        ori_coordinate = [[1, 0, 0, 0] 
                          [0, 1, 0, 0]
                          [0, 0, 1, 0] 
                          [0, 0, 0, 1]];

        new_coordinate(k,:,:) = ori_coordinate * rotate_z(deg2rad(-gravity_direction))* trans(l1, 0, 0) *  rotate_z(rotate_angle(k)) ; 
        new_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [0,0,0,1]';
        
        site4_ = site4 * rotate_z_(deg2rad(gravity_direction));
        site5_ = site5 * rotate_z_(deg2rad(gravity_direction));
        site6_ = site6 * rotate_z_(deg2rad(gravity_direction));

       
        last_point1(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site4 * (kp ^ k), 1]';
        last_point2(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site5 * (kp ^ k), 1]';
        last_point3(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site6 * (kp ^ k), 1]';

        
        active_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site1 * (kp ^ k), 1]';
        passive_point1(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site2 * (kp ^ k), 1]';
        passive_point2(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site3 * (kp ^ k), 1]';
               
        
        passive_tendon_length(k) = dis(site5_, passive_point1(k,:));
        active_tendon_length(k) = dis(site4_, active_point(k,:));
        
    elseif k < section_number % translate and rotate around the new coordinate
        last_coordinate = reshape(new_coordinate(k - 1,:,:), [4,4]);
        new_coordinate(k,:,:) = last_coordinate * trans(l1 * kp ^ (k-1), 0, 0) *  rotate_z(rotate_angle(k));
        new_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [0,0,0,1]'; 
        

        last_point1(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site4 * (kp ^ k), 1]';
        last_point2(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site5 * (kp ^ k), 1]';
        last_point3(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site6 * (kp ^ k), 1]';

        
        active_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site1 * (kp ^ k), 1]';
        passive_point1(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site2 * (kp ^ k), 1]';
        passive_point2(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site3 * (kp ^ k), 1]';
        
       
        passive_tendon_length(k) = dis(last_point2(k - 1,:), passive_point1(k,:));
        active_tendon_length(k) = dis(last_point1(k - 1,:), active_point(k,:));

        
    else
            last_coordinate = reshape(new_coordinate(k - 1,:,:), [4,4]);
            new_coordinate(k,:,:) = last_coordinate * trans(l1 * kp ^ (k - 1), 0, 0);
            new_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [0,0,0,1]';
            
    end
end

t_rest = t_sum;
t_intpu_ini = 37.6085;
t_intpu = t_intpu_ini;
for i = 2 :section_number
    t_intpu = t_intpu + t_intpu_ini * kp ^ (i - 1);
end

passive_tendon_length_sum = sum(passive_tendon_length);
passive_tendon_length_change = passive_tendon_length_sum - t_rest;

active_tendon_length_sum = sum(active_tendon_length);
active_tendon_length_change = t_rest - active_tendon_length_sum;


%% actuating the second cable
% the actuation of the second cable would not affect the first cable. The effect starts from the beginning to start_section
start_section = find(extrem_flag == 1, 1);
new_portion = re_portion(flip(Portion), start_section - 1);

% specify the actuating threshold length from start_section
% consider collision when the cone intersection angle is less than or equal to 150 degrees
threshold_angle = 150;
cone_axis_angle = zeros(1, start_section - 1);
for i = 1 : start_section - 1
    if i == 1
        cone_axis_angle(i) = cone_angle([0,0,0], new_point(i,1:3), new_point(i + 1, 1:3));
    else
        cone_axis_angle(i) = cone_angle(new_point(i - 1, 1:3), new_point(i,1:3), new_point(i + 1, 1:3));
    end
end


new_point_t = zeros(start_section - 1, 3);
optimal_alpha = zeros(start_section - 1, 1);
passive_point1_r = zeros(start_section - 1, 3);
passive_point2_r = zeros(start_section - 1, 3);
active_point_r = zeros(start_section - 1, 3);
most_change_length = zeros(start_section - 1, 1);
minimal = zeros(start_section - 1, 1);
cone_angle_minimal = zeros(start_section - 1, 1);

% A: axis start point
% B: the last joint which is not rotating 
% C: the rotating joint
% D: axis end point
% rotate_vector_to_target_angle_bisection：output the angle where C rotate around AD, which makes angle CAB be a specific angle

for i = start_section - 1 : -1 : 1
    if i > 1
        A = new_point(i,1:3)';
        B = new_point(i - 1, 1:3)';
        C = new_point(i + 1, 1:3)';
        D = last_point1(i - 1, 1:3)';
        [new_point_t(i,:), optimal_alpha(i)] = rotate_vector_to_target_angle_bisection(A, B, C, D, threshold_angle);

        active_point_r(i,:) = rotate_point_around_axis(active_point(i,1:3), A, D, deg2rad(optimal_alpha(i)));
        passive_point1_r(i,:) = rotate_point_around_axis(passive_point1(i,1:3), A, D, deg2rad(optimal_alpha(i)));
        passive_point2_r(i,:) = rotate_point_around_axis(passive_point2(i,1:3), A, D, deg2rad(optimal_alpha(i)));

        cone_angle_minimal(i) = cone_angle(passive_point1_r(i,:), A', last_point2(i - 1, 1:3));

        minimal(i) = dis(passive_point1_r(i,1:3), last_point2(i - 1,1:3));
        most_change_length(i) = passive_tendon_length(i) - minimal(i);

      
    else
        A = new_point(i,1:3)';
        B = [0,0,0]';
        C = new_point(i + 1, 1:3)';
        D = site4_(1:3)';
        [new_point_t(i,:), optimal_alpha(i)] = rotate_vector_to_target_angle_bisection(A, B, C, D, threshold_angle);
        active_point_r(i,:) = rotate_point_around_axis(active_point(i,1:3), A, D, deg2rad(optimal_alpha(i)));
        passive_point1_r(i,:) = rotate_point_around_axis(passive_point1(i,1:3), A, D, deg2rad(optimal_alpha(i)));
        passive_point2_r(i,:) = rotate_point_around_axis(passive_point2(i,1:3), A, D, deg2rad(optimal_alpha(i)));

        cone_angle_minimal(i) = cone_angle(passive_point1_r(i,:), A', site5(1:3));

       
        % calculate the most allowable actuating length to be compared
        minimal(i) = dis(passive_point1_r(i,1:3), site5_(1:3));
        most_change_length(i) = passive_tendon_length(i) - minimal(i);
    end
end

% D2 is the actuation length for the second cable
new_point_r = zeros(start_section - 1, 3);
new_tendon_length2 = passive_tendon_length(1 : start_section - 1);
distribute_length2 = zeros(1, (start_section - 1));
rotate_alpha = zeros(1, (start_section - 1));

for i =(start_section - 1): -1 : 1 
    % current cable length = passive_tendon_length. Distribution length = D2 * new_portion(i)
    % the smallest (most tight) cable length = minimal_length
    % thus, need to compare the largest allowable atcuating length and the distribution length

    distribute_length2(i) = D2 * new_portion(i);
    if distribute_length2(i) < most_change_length(i) 
        new_tendon_length2(i) = passive_tendon_length(i) - distribute_length2(i);
    else % if the distribution length exceed limit, then let the cable length to be the minimal and update the distribution portion
        new_tendon_length2(i) = minimal(i);
        new_portion = re_portion(new_portion, i);
        D2 = D2 - most_change_length(i);
    end
end

% optimal_angle: the exact rotation angle when the second cable actuate most based on the first cable
% cone_angle_minimal: the exact intersection angle when the second cable actuate most based on the first cable
final_angle2 = zeros(1, (start_section - 1));%final_angle2: the exact rotation angle when the second cable is actuated

for k = 1: (start_section - 1)
    final_angle2(k) = tendon2angle2(k, kp, side1, new_tendon_length2(k), hole_point);
    rad2deg(final_angle2(k))
end

active_point_r = zeros(start_section - 1, 3);
passive_point1_r = zeros(start_section - 1, 3);
passive_point2_r = zeros(start_section - 1, 3);
last_point1_r = zeros(start_section - 1, 3);
last_point2_r = zeros(start_section - 1, 3);
last_point3_r = zeros(start_section - 1, 3);

cone_angle2 = zeros(start_section - 1, 1);%cone_angle2: the intersection angle before actuation
for i = 1 : (start_section - 1)
    if i == 1
        A = new_point(i,1:3)';
        B = site5';
        C = passive_point1(i, 1:3)';
        D = site4';
        cone_angle2(i) = cone_angle(C',A',B');
    else  % looking for new rotation angle
        A = new_point(i,1:3)';
        B = last_point2(i - 1, 1:3)';
        C = passive_point1(i, 1:3)';
        D = last_point1(i - 1, 1:3)';
        cone_angle2(i) = cone_angle(C',A',B');
    end
end


%% drawing
figure();
new_point_r(1) = new_point(1);
scatter3(new_point_r(1,1), new_point_r(1,2), new_point_r(1,3));hold on;
plot3([new_point_r(1,1),0], [new_point_r(1,2),0], [new_point_r(1,3),0])

third_line_length = zeros(1, section_number - 1);
second_active_length = zeros(1, section_number - 1);


for i = 1 : (start_section - 1)
    if i == 1
        A = new_point(i,1:3)';
        B = site5_';
        C = passive_point1(i, 1:3)';
        D = site4_';

        % After rotation, the angle between passive_point,A, and last_point2 shall be final_angle2
        % rotate_alpha: the rotation angle where C rotate around axis    
        % the joint new_point(i+1), threading point active_point, passive_point2 also need to rotate

        [passive_point1_r(i,:), rotate_alpha(i)] = rotate_vector_to_target_angle_bisection(A, B, C, D, rad2deg(final_angle2(i)));

        scatter3(site4_(1), site4_(2), site4_(3), 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(site5_(1), site5_(2), site5_(3), 'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(site6_(1), site6_(2), site6_(3), 'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;

        new_point_r(i + 1,:) = rotate_point_around_axis(new_point(i + 1 , 1:3), A, D, deg2rad(rotate_alpha(i)));  
        active_point_r(i,:) = rotate_point_around_axis(active_point(i, 1:3), A, D, deg2rad(rotate_alpha(i)));
        passive_point2_r(i,:) = rotate_point_around_axis(passive_point2(i, 1:3), A, D, deg2rad(rotate_alpha(i)));
        last_point1_r(i,:) = rotate_point_around_axis(last_point1(i , 1:3), A, D, deg2rad(rotate_alpha(i)));  
        last_point2_r(i,:) = rotate_point_around_axis(last_point2(i , 1:3), A, D, deg2rad(rotate_alpha(i)));
        last_point3_r(i,:) = rotate_point_around_axis(last_point3(i , 1:3), A, D, deg2rad(rotate_alpha(i)));

        scatter3(new_point_r(i + 1,1), new_point_r(i + 1,2), new_point_r(i + 1, 3));hold on;

        scatter3(last_point1_r(i,1), last_point1_r(i,2), last_point1_r(i, 3),'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(last_point2_r(i,1), last_point2_r(i,2), last_point2_r(i, 3),'MarkerFaceColor',[0 1 0], 'MarkerEdgeColor',[0 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(last_point3_r(i,1), last_point3_r(i,2), last_point3_r(i, 3),'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor',[0 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
          
        % [1 0 1] rose red   [1 1 0] yellow [0 1 1] sky blue
        % [1 0 0] red [0 1 0] green [0 0 1] blue

        scatter3(active_point_r(i,1), active_point_r(i,2), active_point_r(i, 3),'MarkerFaceColor',[1 0 1], 'MarkerEdgeColor',[1 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(passive_point1_r(i,1), passive_point1_r(i,2), passive_point1_r(i,3),'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor',[1 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(passive_point2_r(i,1), passive_point2_r(i,2), passive_point2_r(i,3),'MarkerFaceColor',[0 1 1], 'MarkerEdgeColor',[0 1 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;

        plot3([active_point_r(i,1),passive_point1_r(i,1),passive_point2_r(i,1),active_point_r(i,1)], ...
            [active_point_r(i,2),passive_point1_r(i,2),passive_point2_r(i,2),active_point_r(i,2)], ...
            [active_point_r(i,3),passive_point1_r(i,3),passive_point2_r(i,3),active_point_r(i,3)]);hold on;
        plot3([new_point_r(1,1),new_point_r(2,1)], [new_point_r(1,2),new_point_r(2,2)], [new_point_r(1,3),new_point_r(2,3)]);

        third_line_length(i) = dis([passive_point2_r(i,1), passive_point2_r(i,2), passive_point2_r(i,3)], ...
            [site6_(i,1), site6_(i,2), site6_(i,3)]);

        second_active_length(i) = dis([passive_point1_r(i,1), passive_point1_r(i,2), passive_point1_r(i,3)], ...
            [site5_(i,1), site5_(i,2), site5_(i,3)]);

        % any point after joint i rotates together
        for k = 2 : (section_number - 1)
            new_point_r(k + 1,:) = rotate_point_around_axis(new_point(k + 1 , 1:3), A, D, deg2rad(rotate_alpha(i)));  
            active_point_r(k,:) = rotate_point_around_axis(active_point(k , 1:3), A, D, deg2rad(rotate_alpha(i)));        
            passive_point1_r(k,:) = rotate_point_around_axis(passive_point1(k , 1:3), A, D, deg2rad(rotate_alpha(i)));
            passive_point2_r(k,:) = rotate_point_around_axis(passive_point2(k , 1:3), A, D, deg2rad(rotate_alpha(i)));

            last_point1_r(k,:) = rotate_point_around_axis(last_point1(k , 1:3), A, D, deg2rad(rotate_alpha(i)));  
            last_point2_r(k,:) = rotate_point_around_axis(last_point2(k , 1:3), A, D, deg2rad(rotate_alpha(i))); 
            last_point3_r(k,:) = rotate_point_around_axis(last_point3(k , 1:3), A, D, deg2rad(rotate_alpha(i))); 

            end

    else % looks for rotation angle
        A = new_point_r(i,1:3)';
        B = last_point2_r(i - 1, 1:3)';
        C = passive_point1_r(i, 1:3)';
        D = last_point1_r(i - 1, 1:3)';
        
        [passive_point1_r(i,:), rotate_alpha(i)] = rotate_vector_to_target_angle_bisection(A, B, C, D, rad2deg(final_angle2(i)));

        new_point_r(i + 1,:) = rotate_point_around_axis(new_point_r(i + 1 , 1:3), A, D, deg2rad(rotate_alpha(i)));  
        active_point_r(i,:) = rotate_point_around_axis(active_point_r(i , 1:3), A, D, deg2rad(rotate_alpha(i)));        
        passive_point2_r(i,:) = rotate_point_around_axis(passive_point2_r(i , 1:3), A, D, deg2rad(rotate_alpha(i)));

        last_point1_r(i,:) = rotate_point_around_axis(last_point1_r(i , 1:3), A, D, deg2rad(rotate_alpha(i)));  
        last_point2_r(i,:) = rotate_point_around_axis(last_point2_r(i , 1:3), A, D, deg2rad(rotate_alpha(i)));
        last_point3_r(i,:) = rotate_point_around_axis(last_point3_r(i , 1:3), A, D, deg2rad(rotate_alpha(i)));

        scatter3(last_point1_r(i,1), last_point1_r(i,2), last_point1_r(i, 3),'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(last_point2_r(i,1), last_point2_r(i,2), last_point2_r(i, 3),'MarkerFaceColor',[0 1 0], 'MarkerEdgeColor',[0 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(last_point3_r(i,1), last_point3_r(i,2), last_point3_r(i, 3),'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor',[0 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;

        scatter3(new_point_r(i + 1,1), new_point_r(i + 1,2), new_point_r(i + 1,3));hold on;
        scatter3(active_point_r(i,1), active_point_r(i,2), active_point_r(i,3),'MarkerFaceColor',[1 0 1], 'MarkerEdgeColor',[1 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(passive_point1_r(i,1), passive_point1_r(i,2), passive_point1_r(i,3),'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor',[1 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(passive_point2_r(i,1), passive_point2_r(i,2), passive_point2_r(i,3),'MarkerFaceColor',[0 1 1], 'MarkerEdgeColor',[0 1 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;

        plot3([active_point_r(i,1),passive_point1_r(i,1),passive_point2_r(i,1),active_point_r(i,1)], ...
            [active_point_r(i,2),passive_point1_r(i,2),passive_point2_r(i,2),active_point_r(i,2)], ...
            [active_point_r(i,3),passive_point1_r(i,3),passive_point2_r(i,3),active_point_r(i,3)]);hold on;

        plot3([new_point_r(i + 1,1),new_point_r(i, 1)], [new_point_r(i + 1,2),new_point_r(i, 2)], [new_point_r(i + 1,3), new_point_r(i, 3)]);
        
        third_line_length(i) = dis([passive_point2_r(i,1), passive_point2_r(i,2), passive_point2_r(i,3)], ...
            [last_point3_r(i - 1,1), last_point3_r(i - 1,2), last_point3_r(i - 1,3)]);

        second_active_length(i) = dis([passive_point1_r(i,1), passive_point1_r(i,2), passive_point1_r(i,3)], ...
            [last_point2_r(i - 1,1), last_point2_r(i - 1,2), last_point2_r(i - 1,3)]);


        for k = (i + 1) : (section_number - 1)
            new_point_r(k + 1,:) = rotate_point_around_axis(new_point_r(k + 1 , 1:3), A, D, deg2rad(rotate_alpha(i)));  
            active_point_r(k,:) = rotate_point_around_axis(active_point_r(k , 1:3), A, D, deg2rad(rotate_alpha(i)));        
            passive_point1_r(k,:) = rotate_point_around_axis(passive_point1_r(k , 1:3), A, D, deg2rad(rotate_alpha(i)));
            passive_point2_r(k,:) = rotate_point_around_axis(passive_point2_r(k , 1:3), A, D, deg2rad(rotate_alpha(i)));
            last_point1_r(k,:) = rotate_point_around_axis(last_point1_r(k , 1:3), A, D, deg2rad(rotate_alpha(i)));  
            last_point2_r(k,:) = rotate_point_around_axis(last_point2_r(k , 1:3), A, D, deg2rad(rotate_alpha(i))); 
            last_point3_r(k,:) = rotate_point_around_axis(last_point3_r(k , 1:3), A, D, deg2rad(rotate_alpha(i))); 

        end
        
        if i == start_section - 1 % regard following sections as rigid body, since the effect of the second cable ends
            for k = (i+1) : (section_number - 1)
                scatter3(new_point_r(k + 1,1), new_point_r(k + 1,2), new_point_r(k + 1,3));hold on;
                scatter3(active_point_r(k,1), active_point_r(k,2), active_point_r(k,3),'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
                scatter3(passive_point1_r(k,1), passive_point1_r(k,2), passive_point1_r(k,3),'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor',[1 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
                scatter3(passive_point2_r(k,1), passive_point2_r(k,2), passive_point2_r(k,3),'MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
                scatter3(last_point1_r(k,1), last_point1_r(k,2), last_point1_r(k, 3),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
                scatter3(last_point2_r(k,1), last_point2_r(k,2), last_point2_r(k, 3),'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
                scatter3(last_point3_r(k,1), last_point3_r(k,2), last_point3_r(k, 3),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
                plot3([active_point_r(k,1),passive_point1_r(k,1),passive_point2_r(k,1),active_point_r(k,1)], ...
            [active_point_r(k,2),passive_point1_r(k,2),passive_point2_r(k,2),active_point_r(k,2)], ...
            [active_point_r(k,3),passive_point1_r(k,3),passive_point2_r(k,3),active_point_r(k,3)]);hold on;
                plot3([new_point_r(k,1),new_point_r(k - 1, 1)], [new_point_r(k,2),new_point_r(k - 1, 2)], [new_point_r(k,3), new_point_r(k - 1, 3)]);
                
                third_line_length(k) = dis([passive_point2_r(k,1), passive_point2_r(k,2), passive_point2_r(k,3)], ...
            [last_point3_r(k - 1,1), last_point3_r(k - 1,2), last_point3_r(k - 1,3)]);

                second_active_length(i) = dis([passive_point1_r(i,1), passive_point1_r(i,2), passive_point1_r(i,3)], ...
            [last_point2_r(k - 1,1), last_point2_r(k - 1,2), last_point2_r(k - 1,3)]);
            end
        end
    end
end

xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

%% mujoco simulation
t_intpu_ini = 37.6085;
t_intpu = t_intpu_ini;

for i = 2 :section_number
    t_intpu = t_intpu + t_intpu_ini * kp ^ (i - 1);
end

third_line_sum = sum(third_line_length);
third_line_change = t_sum - third_line_sum;

second_active_sum = sum(second_active_length);
second_active_change = t_sum - second_active_sum;


mujoco_active1 = t_intpu + active_tendon_length_sum
mujoco_active2 = t_intpu + passive_tendon_length_sum
mujoco_passive = t_intpu + third_line_sum

%% functions

function angle = cone_angle(A,B,C)
    % vector BA and BC
    BA = A - B; 
    BC = C - B; 
    
    dot_product = dot(BA, BC); 
    norm_BA = norm(BA); 
    norm_BC = norm(BC); 
    
    theta_rad = acos(dot_product / (norm_BA * norm_BC));
    
    angle = rad2deg(theta_rad);
end

function tendon = angle2tendon(a, b, angle)
    tendon = sqrt(a ^ 2 + b ^ 2 - 2 * a * b * cos(angle));
end

function T_z = rotate_z_(angle)
    T_z = [
    [cos(angle),  -sin(angle), 0]
    [sin(angle), cos(angle), 0]
    [0, 0, 1]
    ];
end

function [AC_rotated, optimal_alpha] = rotate_vector_to_target_angle_bisection(A, B, C, D, target_angle)
    % A, B, C: position
    % D: the direction vector of the rotation axis
    % target_angle: target angle in degree

    AB = B - A;
    AC = C - A;

    initial_angle = acos(dot(AB, AC) / (norm(AB) * norm(AC))) * (180 / pi);

    % if initial angle equals to target angle, then return
    if abs(initial_angle - target_angle) < 1e-6
        AC_rotated = AC;
        optimal_alpha = 0;
        return;
    end

    % initialize the bound for the dichotomy
    lower_bound = - pi / 4;
    upper_bound = 0; 

    tolerance = 1e-6; 

    % dichotomy search for optimal angle
    while (upper_bound - lower_bound) > tolerance
        % select mid point
        alpha = (lower_bound + upper_bound) / 2;
        AC_test = rotate_point_around_axis(C, A, D, alpha);
        current_angle = acos(dot(AB, AC_test - A) / (norm(AB) * norm(AC_test - A))) * (180 / pi);

        % adjust bound
        if current_angle > target_angle
            upper_bound = alpha;
        else
            lower_bound = alpha;
        end
    end

    optimal_alpha = (lower_bound + upper_bound) / 2;
    AC_rotated = rotate_point_around_axis(C, A, D, optimal_alpha);
    optimal_alpha = rad2deg(optimal_alpha);
end

function new_portion = re_portion(portion, index)
    new_portion = zeros(1, length(index));
    new_sum = sum(portion(1 : index));
    for i = 1 : index
        new_portion(i) = portion(i) / new_sum;
    end
end

function p_rotated = rotate_point_around_axis(p_to_rotate, axis_start, axis_end, rotation_angle)
    % Convert points to vectors
    p_to_rotate = p_to_rotate(:);
    axis_start = axis_start(:);
    axis_end = axis_end(:);

    % Find the direction vector of the axis
    axis_dir = axis_end - axis_start;

    % Normalize the direction vector
    axis_normal = axis_dir / norm(axis_dir);

    % Translate point P so that the axis passes through the origin
    P_prime = p_to_rotate - axis_start;

    % Create the skew-symmetric matrix K
    K = [0, -axis_normal(3), axis_normal(2);
         axis_normal(3), 0, -axis_normal(1);
         -axis_normal(2), axis_normal(1), 0];

    % Calculate the rotation matrix R using Rodrigues' rotation formula
    I = eye(3); % Identity matrix
    R = I + sin(rotation_angle) * K + (1 - cos(rotation_angle)) * (K * K);

    % Rotate the translated point
    P_rotated = R * P_prime;

    % Translate back to the original position
    p_rotated = P_rotated + axis_start;
end

% 3D translation matrix
function T_t = trans(x, y, z)
    T_t = [[1, 0, 0, x]
           [0, 1, 0, y]
           [0, 0, 1, z]
           [0, 0, 0, 1]];
end


% rotation around z axis
function T_z = rotate_z(angle)
    T_z = [
    [cos(angle),  -sin(angle), 0, 0]
    [sin(angle), cos(angle), 0, 0]
    [0, 0, 1, 0]
    [0, 0, 0, 1]];
end

function angle = tendon2angle(a, b, tendon_length)
    if tendon_length == 0
        angle = 0;
    else
        cos_angle = (a ^ 2 + b ^ 2 - tendon_length ^ 2) / (2 * a * b);
        angle = acos(cos_angle);
    end
end

function angle = tendon2angle2(m, kp, side1, tendon_length, hole_point)
    if tendon_length == 0
        angle = 0;
    else
        side = side1 * kp ^ (m - 1);
        a = hole_point * side;
        b = hole_point * side;
        cos_angle = (a ^ 2 + b ^ 2 - tendon_length ^ 2) / (2 * a * b);
        angle = acos(cos_angle);
    end
end

function distance = dis(p1, p2)
    distance = sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2 + (p1(3) - p2(3))^2);
end

end