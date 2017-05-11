function [ marker_estimate ] = linearize_C( state_estimate, visible_markers)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global n_m markers_body

marker_estimate = zeros(length(visible_markers),1);
valid_ind = 1:3;

for i = 1:n_m
    ind = (1:3) + (i-1)*3;
    if(any(visible_markers == ind(1)))
        marker_estimate(valid_ind) = state_estimate(1:3) + quat2rotm(state_estimate(7:10)')*markers_body(i,:)';
        valid_ind = valid_ind + 3;
    end
end

end

