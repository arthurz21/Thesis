function [qdur_avF] = qdur_avF(ecg, QRST_pts)

    % Extract fiducial points
    Q_start_pt = QRST_pts(1);
    R_peak_estimate = QRST_pts(2);
    S_end_pt = QRST_pts(3);

    % Find the R peak magnitude within the range Q_start to S_end
    [R_peak_mag, relative_index] = max(ecg.avF(Q_start_pt:S_end_pt));

    % Compute the global index of the R peak
    R_peak_index = Q_start_pt + relative_index - 1;
    Q_end = find_R_wave_start_avF(ecg, R_peak_index, Q_start_pt, R_peak_mag);
    
    qdur_avF = Q_end - Q_start_pt;
end


function Q_end = find_R_wave_start_avF(ecg, R_peak_index, Q_start_pt, R_peak_mag)
    
    % Define lead names (standard 12-lead ECG)
    
    % Get fiducial points
    Q_start = Q_start_pt;
    R_peak = R_peak_index;
        
    % Initialize R_wave_start as NaN
    Q_end = NaN;

    % Extract the data for the current lead
    lead_data = ecg.avF;

    % Calculate 20% threshold of R peak magnitude
    threshold = 0.2 * R_peak_mag;

    % Initialize variables for both methods
    zero_cross_point = NaN;
    threshold_point = NaN;

     %% Find R_wave_start by searching backward from R_peak
     for j = R_peak:-1:Q_start+1
         % Check for sign change between consecutive points
         if sign(lead_data(j)) ~= sign(lead_data(j - 1))
             zero_cross_point = j; % Set R_wave_start to the first sign change
             break;
         end
     end

    % Search for threshold crossing
    for j = R_peak:-1:Q_start+1
        if abs(lead_data(j)) <= threshold
            threshold_point = j;
            break;
        end
    end

    % Decision logic for combining methods
    if ~isnan(zero_cross_point) && ~isnan(threshold_point)
        % Use the earlier point (closer to QRS onset)
        Q_end = min(zero_cross_point, threshold_point);
    elseif ~isnan(zero_cross_point)
        Q_end = zero_cross_point;
    elseif ~isnan(threshold_point)
        Q_end = threshold_point;
    end

end
