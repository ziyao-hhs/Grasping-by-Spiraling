% gravity_angle: right horizontal 0°, vertial 0°, left horizontal 180°

%% dorsal single bend
% input: cable actuation (mm), gravity angle (°)
% output: the length change of the second and third passive cable & the predicted joint position

[passive_tendon_length_change_dorsal, joint_dorsal] = dorsal_single_bend_fun(50, 90);

%% ventral two cable bend
% input: cable actuation (mm), gravity angle (°)
% output: the length change of the second and third passive cable & the predicted joint position

% [passive_tendon_length_change_ventral, joint_ventral] = ventral_two_cable_bend_fun(100, 90);

%% twisting
% input: first cable actuation (mm), second cable actuation (mm), gravity angle (°)
% output: the length change of the third passive cable & the predicted joint position

% [passive_tendon_length_change_twsit, joint_twist] = twist_gravity_fun(60, 30, 90);
