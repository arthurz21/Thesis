function [global_first_inflection_amplitude, fpTinfl1Axis, fpTinfl1Mag]= calculate_global_first_inflection(ecg, T_start, T_end)
    % Inputs:
    % - ecg: A struct with fields 'X', 'Y', and 'Z' representing the ECG signals
    % - T_points: A 3x2 matrix where each row corresponds to [T_start, T_end] for each lead

    % Initialize inflection points for each lead as NaN
    XT_first_inflection = NaN;
    YT_first_inflection = NaN;
    ZT_first_inflection = NaN;

    % Function to calculate the second derivative
    second_derivative = @(data) diff(data, 2);

    % Loop over each lead to find the first inflection point between T_start and T_end
    lead_names = {'X', 'Y', 'Z'};
    inflection_points = [XT_first_inflection, YT_first_inflection, ZT_first_inflection];
    
    
    for i = 1:length(lead_names)
        lead = lead_names{i};
        lead_data = ecg.(lead);
        

        % Calculate the second derivative
        sec_derivative = second_derivative(lead_data);

        % Find the first index between T_start and T_end where the second derivative is near zero
        for j = T_start:T_end-1
            if abs(sec_derivative(j)) < 1e-3  % Threshold for "nearing zero"
                inflection_points(i) = j;
                break;
            end
        end
    end
    
    XT_first_inflection = inflection_points(1);
    YT_first_inflection = inflection_points(2);
    ZT_first_inflection = inflection_points(3);

    % Check if each inflection point is valid and get the amplitude from the original ECG data
    if ~isnan(XT_first_inflection)
        XT_amplitude = ecg.X(XT_first_inflection);
    end
    if ~isnan(YT_first_inflection)
        YT_amplitude = ecg.Y(YT_first_inflection);
    end
    if ~isnan(ZT_first_inflection)
        ZT_amplitude = ecg.Z(ZT_first_inflection);
    end

    % Compute the global amplitude using the valid amplitudes
    global_first_inflection_amplitude = sqrt(...
        (XT_amplitude * ~isnan(XT_amplitude))^2 + ...
        (YT_amplitude * ~isnan(YT_amplitude))^2 + ...
        (ZT_amplitude * ~isnan(ZT_amplitude))^2 ...
    );
    
    fpTinfl1Mag = sqrt(XT_amplitude*XT_amplitude + YT_amplitude*YT_amplitude)*1000;
    fpTinfl1Axis = atan2(YT_first_inflection, XT_first_inflection);
    fpTinfl1Axis = rad2deg(fpTinfl1Axis);

    % If all amplitudes are NaN, set global amplitude to NaN
    if isnan(XT_amplitude) && isnan(YT_amplitude) && isnan(ZT_amplitude)
        global_first_inflection_amplitude = NaN;
    end
    

end
