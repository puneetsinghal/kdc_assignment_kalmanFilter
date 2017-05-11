function [ marker_error, C ] = predict_markers( state_estimate, index, markersn )
%
%
global n_m markers_body

marker_error = zeros(n_m*3,1);
marker_estimate = zeros(n_m*3,1);

for i = 1:n_m
  ind = (1:3) + (i-1)*3;
  marker_estimate(ind) = state_estimate(1:3) + quat2rotm(state_estimate(7:10)')*markers_body(i,:)';
  marker_error(ind) = markersn(index, ind)' - marker_estimate(ind);
end

%% Numerical Approach to finding Jacobian Matrix
epsilon = 1e-6;
C = zeros(n_m*3,length(state_estimate));

for i = 1:length(state_estimate)
    state_new = state_estimate;
    state_new(i) = state_new(i) + epsilon;
    new_marker_estimate = linearize_C(state_new);
    
    state_old = state_estimate;
    state_old(i) = state_old(i) - epsilon;
    old_marker_estimate = linearize_C(state_old);
    C(:,i) = (new_marker_estimate - old_marker_estimate)/(2*epsilon);
end

end

