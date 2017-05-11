function [ state_estimate ] = linearize_A( state_estimate, Iworld )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global T

q0 = state_estimate(7);
q1 = state_estimate(8);
q2 = state_estimate(9);
q3 = state_estimate(10);
wx = state_estimate(11);
wy = state_estimate(12);
wz = state_estimate(13);

I_body_x = 1.571428571428571;
I_body_y = 5.362637362637362;
I_body_z = 7.065934065934067;
Ibody = diag([I_body_x, I_body_y, I_body_z]);

% Calculation for Angular Acceleration
skew_omega = hat(state_estimate(11:13));

Rq = quat2rotm(state_estimate(7:10)');
[U,~,V] = svd(Rq);
Rq = U*V';

Ra = [1 -wz*T wy*T; wz*T 1 -wx*T; -wy*T wx*T 1];
[U,~,V] = svd(Ra);
Ra = U*V';

Iworld = Ra*Rq*Ibody*Rq'*Ra';

wd = -Iworld\skew_omega*Iworld*[wx; wy; wz;];

% UPDATE QUATERNION
q0_new = q0 + T*1/(2*sqrt(q0^2 + q1^2 + q2^2 + q3^2))*(-wx*q1 - wy*q2 - wz*q3);
q1_new = q1 + T*1/(2*sqrt(q0^2 + q1^2 + q2^2 + q3^2))*(wx*q0 + wy*q3 - wz*q2);
q2_new = q2 + T*1/(2*sqrt(q0^2 + q1^2 + q2^2 + q3^2))*(wy*q0 + wz*q1 - wx*q3);
q3_new = q3 + T*1/(2*sqrt(q0^2 + q1^2 + q2^2 + q3^2))*(wz*q0 + wx*q2 - wy*q1);

quaternion = [q0_new; q1_new; q2_new; q3_new];
quaternion = quaternion/norm(quaternion);


%% Update State Estimate

% COM POSITION
for i = 1:3
 state_estimate(i) = state_estimate(i) + T*state_estimate(i+3);
end

% LINEAR VELOCITY REMAINS SAME UNDER ZERO FORCE
state_estimate(4:6) = state_estimate(4:6);

% Quaternion
state_estimate(7:10) = quaternion;

% Angular Velocity
state_estimate(11:13) = state_estimate(11:13) + T*wd;

end

