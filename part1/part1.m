clear; clc;
addpath(genpath(fullfile(pwd,'p1n')));

figure;

for filesIndex = 0:9
    clc;
    fileName = strcat('p1n0', int2str(filesIndex));
    markersn = load(fullfile(pwd, 'p1n', fileName));
    
    % load marker data file
%     load markersn_ass4.mat
    
    global T
    % time step (length of time between measurements)
    T = 0.01;
    
    % number of filter states
    global n_x
    n_x = 13;
    
    % number of markers
    global n_m
    n_m = 8;
    
    % number of samples
    n_t = length(markersn);
    
    % marker locations in body coordinates
    global markers_body
    
    markers_body = [
        -1.5   -3.0   -4.5
        -1.5   -3.0    1.5
        -1.5    1.0   -4.5
        -1.5    1.0    1.5
        0.5   -3.0   -4.5
        0.5   -3.0    1.5
        0.5    1.0   -4.5
        0.5    1.0    1.5
        ];
    
    % markers_body = (rotz(0)*markers_body')';
    
    %% set up data arrays
    quaternions_before_measurement( n_t, 4 ) = 0;
    states_before_measurement( n_t, n_x ) = 0;
    quaternions_before_process_update( n_t, 4 ) = 0;
    states_before_process_update( n_t, n_x ) = 0;
    marker_errors( n_t, n_m*3 ) = 0;
    kalman_gains( n_t, 3 ) = 0;
    variance_estimates_before_measurement( n_t, n_x ) = 0;
    variance_estimates_before_process_update( n_t, n_x ) = 0;
    
    %% Initial estimates. Since EKF, need to get nonlinear part to be close
    [ quaternion, state_estimate ] = hack_init_estimates( markersn );
    % state_estimate = [-11 0 0 0 0 0]';
    % state_estimate(4:6) = [0.05, 0.02, 0.01];
    
    % Initial condition noise model
    variance_estimate = eye( n_x );
    
    % postions
    variance_estimate(1, 1) = 1;
    variance_estimate(2, 2) = 1;
    variance_estimate(3, 3) = 1;
    
    % velocities
    variance_estimate(4, 4) = 1;
    variance_estimate(5, 5) = 1;
    variance_estimate(6, 6) = 1;
    
    variance_estimate(7, 7) = 1;
    variance_estimate(8, 8) = 1;
    variance_estimate(9, 9) = 1;
    variance_estimate(10, 10) = 1;
    
    variance_estimate(11, 11) = 1;
    variance_estimate(12, 12) = 1;
    variance_estimate(13, 13) = 1;
    % variance_estimate = variance_estimate;
    
    Q = diag( ones( n_x, 1 ) );
    Q(1,1) = 1*1e-10;
    Q(2,2) = 1*1e-10;
    Q(3,3) = 1*1e-10;
    Q(4,4) = 1e-10;
    Q(5,5) = 1e-10;
    Q(6,6) = 1e-10;
    Q(7,7) = 1*1e-10;
    Q(8,8) = 1*1e-10;
    Q(9,9) = 1*1e-10;
    Q(10,10) = 1*1e-10;
    Q(11,11) = 1e-8;
    Q(12,12) = 1e-8;
    Q(13,13) = 1e-8;
    
    
    % Measurement noise model: marker variance
    R = 0.1*diag(ones(n_m*3, 1));
    
    %%
    for index = 1:length(markersn)
        
        % Relinearize
        %  [ quaternion, state_estimate ] = relinearize( quaternion, state_estimate );
        
        % Compute measurement error
        [ marker_error, C ] = predict_markers( state_estimate, ...
            index, markersn );
        
        % Store values
        quaternions_before_measurement( index, : ) = quaternion;
        states_before_measurement( index, : ) = state_estimate;
        marker_errors( index, : ) = marker_error;
        for i = 1:n_x
            variance_estimates_before_measurement( index, i ) = ...
                variance_estimate( i, i );
        end
        
        % Measurement update
        % should not use inv(), solve set of equations instead.
        kalman_gain = variance_estimate*C'/(C*variance_estimate*C' + R);
        
        state_estimate = state_estimate + kalman_gain*marker_error;
        variance_estimate = (eye(n_x) - kalman_gain*C)*variance_estimate;
        
        % Relinearize
        %  [ quaternion, state_estimate ] = relinearize( quaternion, state_estimate );
        
        % Store values
        quaternions_before_process_update( index, : ) = quaternion;
        states_before_process_update( index, : ) = state_estimate;
        for i = 1:n_x
            variance_estimates_before_process_update( index, i ) = ...
                variance_estimate( i, i );
        end
        
        % Process update.
        [ state_estimate, A ] = do_dynamics( state_estimate );
        variance_estimate = A*variance_estimate*A' + Q;
        
    end
    
    %% Saving file
    fid = fopen(fullfile(pwd, 'results', strcat('p1a0', int2str(filesIndex), '.txt')), 'wt');
    fprintf(fid, 'x\ty\tz\txd\tyd\tzd\tq0\tq1\tq2\tq3\twx\twy\twz\n');
    for i = 1:n_t
        fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',...
            states_before_measurement(i,1), states_before_measurement(i,2), states_before_measurement(i,3),...
            states_before_measurement(i,4), states_before_measurement(i,5), states_before_measurement(i,6),...
            states_before_measurement(i,7), states_before_measurement(i,8), states_before_measurement(i,9), states_before_measurement(i,10),...
            states_before_measurement(i,11), states_before_measurement(i,12), states_before_measurement(i,13));
    end
    fclose(fid);
    
    %% Plotting
    truth_1 = load('problem_2_0.dat');
    truth_2 = load('problem_2_1.dat');
    close all
    time = linspace(0,n_t*T,n_t);
    figure;
    
    subplot(2,2,1);
    hold on
    plot(time, states_before_process_update(:,11),'r',...
        time, states_before_process_update(:,12),'g',...
        time, states_before_process_update(:,13),'b');
%     plot(time, truth_2(1:n_t,1),'r',...
%         time, truth_2(1:n_t,2),'g',...
%         time, truth_2(1:n_t,3),'b');
    title('Angular Velocity')
    xlabel('time (s)')
    ylabel('rad/s')
    legend('x-Angular Velocity','y-Angular Velocity','z-Angular Velocity','Location','best');
    
    
    subplot(2,2,2);
    hold on
    plot(time, states_before_process_update(:,7),'r',...
        time, states_before_process_update(:,8),'g',...
        time, states_before_process_update(:,9),'b',...
        time, states_before_process_update(:,10),'k');
%     plot(time, truth_1(1:n_t,4),'r',...
%         time, truth_1(1:n_t,5),'g',...
%         time, truth_1(1:n_t,6),'b',...
%         time, truth_1(1:n_t,7),'k');
    title('Quaternions')
    xlabel('time (s)')
    legend('q-scalar','q1','q2','q3','Location','best');
    
    subplot(2,2,3)
    hold on
    plot(time, states_before_process_update(:,1),'r',...
        time, states_before_process_update(:,2),'g',...
        time, states_before_process_update(:,3),'b');
%     plot(time, truth_1(1:n_t,1),'r',...
%         time, truth_1(1:n_t,2),'g',...
%         time, truth_1(1:n_t,3),'b');
    title('Center of mass Position')
    xlabel('time (s)')
    ylabel('meters')
    legend('COM Position-X','COM Position-Y','COM Position-Z','Location','best');
    
    subplot(2,2,4)
    hold on
    plot(time, states_before_process_update(:,4),'r',...
        time, states_before_process_update(:,5),'g',...
        time, states_before_process_update(:,6),'b');
%     plot(time, 0.02*ones(n_t,1),'r--',...
%         time, -0.05*ones(n_t,1),'g--',...
%         time, 0.01*ones(n_t,1),'b--');
    title('Center of mass Velocity')
    xlabel('time (s)')
    ylabel('m/s')
    legend('COM Velocity-X','COM Velocity-Y','COM Velocity-Z','Location','best');
    
    saveas(gcf, fullfile(pwd, 'results', strcat('p1a0', int2str(filesIndex), '.png')));
end