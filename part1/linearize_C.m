function [ marker_estimate ] = linearize_C( state_estimate )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global n_m markers_body

marker_error = zeros(n_m*3,1);
marker_estimate = zeros(n_m*3,1);

for i = 1:n_m
  ind = (1:3) + (i-1)*3;
  marker_estimate(ind) = (quat2rotm(state_estimate(7:10)')*markers_body(i,:)') + state_estimate(1:3);
end

end

