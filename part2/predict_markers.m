function [ marker_error, C ] = predict_markers( state_estimate, index, markersn, visible_markers )
%
%
global n_m markers_body

marker_error = zeros(length(visible_markers),1);
marker_estimate = zeros(length(visible_markers),1);
valid_ind = 1:3;

for i = 1:n_m
    ind = (1:3) + (i-1)*3;
    if(any(visible_markers == ind(1)))
        marker_estimate(valid_ind) = state_estimate(1:3) + quat2rotm(state_estimate(7:10)')*markers_body(i,:)';
        marker_error(valid_ind) = markersn(index, ind)' - marker_estimate(valid_ind);
        valid_ind = valid_ind + 3;
    end
end

%% Numerical Approach to finding Jacobian Matrix
epsilon = 1e-6;
C = zeros(length(visible_markers),length(state_estimate));

for i = 1:length(state_estimate)
    state_new = state_estimate;
    state_new(i) = state_new(i) + epsilon;
    new_marker_estimate = linearize_C(state_new, visible_markers);
    
    state_old = state_estimate;
    state_old(i) = state_old(i) - epsilon;
    old_marker_estimate = linearize_C(state_old, visible_markers);
    C(:,i) = (new_marker_estimate - old_marker_estimate)/(2*epsilon);
end

end

