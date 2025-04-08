
function [R_wave_limits, J_wave_limits, ST_segment_limits] = index_finder_vcg(medianvcg1, fidpts)
%cut out ST segments. we will also be checking for Jwave_here
%we already know S
R_wave_limits = find_R_wave_start_end_with_sign_check(medianvcg1, fidpts);
J_wave_limits = find_J_wave_and_ST_start_end(medianvcg1, fidpts);
ST_segment_limits = find_ST_seg(medianvcg1, fidpts);

end


function R_wave_limits = find_R_wave_start_end_with_sign_check(ecg, fidpts)
    lead_names = {'X', 'Y', 'Z'};
    
    % Get fiducial points
    QRS_start = fidpts(1);
    R_peak = fidpts(2);
    S_end = fidpts(3);
    
    % Initialize arrays to store R wave start and end values for averaging
    R_wave_start_values = NaN(1, length(lead_names)); % NaN array for start averaging
    R_wave_end_values = NaN(1, length(lead_names));   % NaN array for end averaging

    % Initialize output structure
    R_wave_limits = struct('avg_start', NaN, 'avg_end', NaN);

    % Loop over each lead to calculate R_wave_start and R_wave_end
    for i = 1:length(lead_names)
        lead = lead_names{i};
        
        % Initialize R_wave_start and R_wave_end as NaN
        R_wave_start = NaN;
        R_wave_end = NaN;

        % Extract the data for the current lead
        lead_data = ecg.(lead);

        %% Find R_wave_start by searching backward from R_peak
        for j = R_peak:-1:QRS_start+1
            % Check for sign change between consecutive points
            if sign(lead_data(j)) ~= sign(lead_data(j - 1))
                R_wave_start = j; % Set R_wave_start to the first sign change
                break;
            end
        end

        %% Find R_wave_end by searching forward from R_peak
        for j = R_peak:S_end-1
            % Check for sign change between consecutive points
            if sign(lead_data(j)) ~= sign(lead_data(j + 1))
                R_wave_end = j; % Set R_wave_end to the point after the sign change
                break;
            end
        end

        % Store the results in arrays for averaging
        R_wave_start_values(i) = R_wave_start;
        R_wave_end_values(i) = R_wave_end;
    end

    %% Calculate Averaged Start and End Points Across All Leads
    % Take the mean of all detected start and end points, omitting NaN values
    avg_start = round(mean(R_wave_start_values, 'omitnan'));
    avg_end = round(mean(R_wave_end_values, 'omitnan'));

    % Store the averaged results in the output structure
    R_wave_limits.avg_start = avg_start;
    R_wave_limits.avg_end = avg_end;

    % Display the results
    %disp('R_wave_start and R_wave_end for each lead:');
    %disp(R_wave_limits);
    %fprintf('Averaged R_wave_start: %d, Averaged R_wave_end: %d\n', avg_start, avg_end);
end


function J_wave_limits = find_J_wave_and_ST_start_end(ecg, fidpts)
    % This function detects the J wave using a larger moving average window to reduce noise
    % Inputs: 
    %   ecg - structure containing the ECG signal with a 'VM' lead for analysis
    %   fidpts - array with key fiducial points [QRS_start, R_peak, S_end, T_peak]
    
    % Parameters
    lead_name = 'VM';            % Use VM lead for J wave detection
    moving_avg_window = 10;      % Window size for the moving average
    derivative_window = 10;      % Window size to detect consistent trends

    % Get fiducial points
    S_end = fidpts(3);           % S wave end point, marking the start of the J wave check
    T_peak = fidpts(4);          % T wave peak, marking the end of the J wave check region

    % Extract the relevant segment of the VM lead
    vm_data = ecg.(lead_name);
    segment_data = vm_data(S_end:T_peak);

    % Calculate derivatives and apply a larger moving average
    derivatives = diff(segment_data);                 % First derivative (rate of change)
    smoothed_derivatives = movmean(derivatives, moving_avg_window);  % Smooth the derivative
    
    % Initialize flags and default values for J wave detection
    J_wave_start = S_end;
    J_wave_end = NaN;  % Assume J wave is not detected unless criteria are met
    marker_rising = false;  % Rising trend marker
    marker_falling = false; % Falling trend marker

    % Loop through the smoothed derivative values to find the J wave
    for i = 1:(length(smoothed_derivatives) - derivative_window + 1)
        % Define a larger window to detect consistent trends
        window = smoothed_derivatives(i:i + derivative_window - 1);

        % Check if all values in the window are consistently positive (rising trend)
        if all(window > 0)
            marker_rising = true;
        elseif all(window < 0) && marker_rising  % Falling trend after a rising trend
            marker_falling = true;
            J_wave_end = S_end + i;  % Mark the end of the J wave
            break;
        end
    end

    % Store the detected J wave start and end points
    J_wave_limits.start = J_wave_start;
    J_wave_limits.end = J_wave_end;

    % Display results
    %disp('J wave limits:');
    %disp(J_wave_limits);
end


function ST_segment_limits = find_ST_seg(ecg, fidpts)
    % Using VM lead to detect ST_end (start of T wave)
    lead_name = 'VM';

    % Get fiducial points
    S_end = fidpts(3);  % ST segment starts at S_end
    T_peak = fidpts(4); % T wave peak
    
    % Parameters
    window_size = 10;              % Window size to check consecutive derivatives
    moving_avg_window = 10;       % Window size for moving average
    vm_data = ecg.(lead_name);
    
    % Compute derivatives
    derivatives = diff(vm_data(S_end:T_peak));
    
    % Apply moving average to derivatives
    smoothed_derivatives = movmean(derivatives, moving_avg_window);
    
    % Initialize tracking variables
    ST_end = NaN; % This will mark the start of T wave
    in_rising_trend = false;

    % Loop through the smoothed derivatives to detect consistent rising trend
    for i = 1:(length(smoothed_derivatives) - window_size + 1)
        % Define the current window of smoothed derivatives
        window = smoothed_derivatives(i:i + window_size - 1);
        
        % Check if we have a consistent rising trend in this window
        is_rising = all(window > 0); % Rising with positive slope
        
        if is_rising
            % If this is the start of a rising trend, mark as potential ST_end (T wave start)
            if ~in_rising_trend
                ST_end = S_end + i - 1; % Calculate ST_end index
                in_rising_trend = true;
            end
        elseif in_rising_trend
            % If we previously had a rising trend, but now it stops,
            % we assume we've identified the ST_end and can exit
            break;
        end
    end
    
    % Store the detected ST_end
    ST_segment_limits.ST_end = ST_end; % End of ST segment (start of T wave)
    ST_segment_limits.ST_start = S_end; % Start of ST segment
    
    % Display results
    disp('ST segment start (S_end):');
    disp(ST_segment_limits.ST_start);
    disp('ST segment end (T wave start detected at):');
    disp(ST_segment_limits.ST_end);
end