function thetadeg = solve_angle(Init_Theta)

%% constant parameters
E = 16.9e6;
Kc = 26.749 / 47.94;
Kp = 0.9196;
g = 9.81;
n = 17;
Theta = sym('Theta', [1 n]);

%% initialize the zero and first section parameters of the trunk
a0 = 52.1271 * 1e-3;
R1 = 3.825 * 1e-3;
L1 = 2.10 * 1e-3;

%% loop for other section parameters
a = zeros(1, n);
R = zeros(1, n);
L = zeros(1, n);

for j = 1 : n
    a(j) = a0 * Kp ^ j;
    R(j) = R1 * Kp ^ (j - 1);
    L(j) = L1 * Kp ^ (j - 1);
end    

m_gripper = 3.2;

% the weight for each section
m = [26, 21.3, 20.4, 14.2, 11.3, 9.7, 7.7, 9.7, 8.5, 7, 5.6, 4.7, 3.5, 3, 2.3, 2, 1.4 + m_gripper] * 1e-3;


%% list force equillibrium formulas
eqn = sym(zeros(1, n));

for i = 1 : n
    if i == n
        eqn(i) = E * I_Moment(R(i)) * Theta(i) / L(i) - m(i) * g * a(i) * Kc * cos(Init_Theta + sum(Theta(1:i)));
    else 
        eqn(i) = E * I_Moment(R(i)) * Theta(i) / L(i) - E * I_Moment(R(i + 1)) * Theta(i + 1) / L(i + 1) - (m(i) * Kc + sum(m(i + 1:n))) * g * a(i) * cos(Init_Theta + sum(Theta(1:i)));
    end
end

%% solve eqn
sol = vpasolve(eqn, Theta);

soltheta = zeros(1,n + 1);
% soltheta(1) = Init_Theta;
soltheta1 = zeros(1,n);

A = zeros(1, n);
A(1) = 52.1271;

for j = 2 : n + 1
    A(j) = A(1) * Kp ^ (j - 1);
end 

P = zeros(2, n + 1);   
% P(:, 1) = [A(1) * cos(Init_Theta); - A(1) * sin(Init_Theta)];


%% output angles
for i = 1 : n
    soltheta1(i) = eval('sol.Theta' + string(i));
end

thetadeg = rad2deg(soltheta1);


function I = I_Moment(R)
    I = pi / 4 * R^4;
end

end











