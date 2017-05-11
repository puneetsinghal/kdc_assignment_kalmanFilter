function [ quaternion, state_estimate ] = hack_init_estimates( markersn )
%
%
global n_m markers_body

marker_rotated = (rotz(-pi/2)*markers_body')';
% use linear regression to estimate COM location
A = zeros( n_m*3, 3 );
b = zeros(n_m*3,1);

% load up matrices: A*p = b
for i = 1:n_m
 for j = 1:3
  row = (i - 1) * 3;
  A(row + j, j) = 1;
  b(row + j) = markersn( 1, row + j ) - marker_rotated(i, j);
 end
end

% actually do regression
p = A\b;
% check errors (residuals)
A*p - b;

quaternion = rand(4,1);
quaternion = quaternion/norm(quaternion);
quaternion = [1 0 0 0]';
state_estimate = [ p(1) p(2) p(3) 0 0 0 quaternion' 0 0 0]';
% state_estimate(1:3) = [11; 0; 4];
