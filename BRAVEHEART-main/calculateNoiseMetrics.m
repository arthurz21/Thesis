function [rms_min_all, rms_var_all] = calculateNoiseMetrics(ecg_raw_cell, filtered_ecg_cell)
    % Calculate noise metrics by comparing raw and filtered ECGs
    %
    % Input:
    %   ecg_raw_cell - Cell array of raw ECGs (ECG12 objects)
    %   filtered_ecg_cell - Cell array of filtered ECGs (ECG12 objects)
    %
    % Output:
    %   rms_min_all - RMS of last 25 samples difference (flat spot noise)
    %   rms_var_all - RMS of entire signal difference (total noise variance)
    
    num_patients = length(ecg_raw_cell);
    rms_min_all = zeros(num_patients, 1);
    rms_var_all = zeros(num_patients, 1);
    
    for patient = 1:num_patients
        % Get current patient's raw and filtered ECGs
        raw = ecg_raw_cell{patient};
        filtered = filtered_ecg_cell{patient};
        
        % Get the 8 independent leads
        raw_leads = {raw.I, raw.II, raw.V1, raw.V2, raw.V3, raw.V4, raw.V5, raw.V6};
        filt_leads = {filtered.I, filtered.II, filtered.V1, filtered.V2, filtered.V3, filtered.V4, filtered.V5, filtered.V6};
        
        % For flat spot noise (last 25 samples)
        flat_rms = zeros(8, 1);
        for i = 1:8
            raw_signal = raw_leads{i}*1000;
            filt_signal = filt_leads{i}*1000;
            
            % Get last 15 samples difference
            noise_diff = raw_signal(end-24:end) - 0; %filt_signal(end-24:end);
            flat_rms(i) = sqrt(mean(noise_diff.^2));
        end
        
        % For total noise variance (all samples)
        noise_var = zeros(8, 1);
        for i = 1:8
            raw_signal = raw_leads{i}*1000;
            filt_signal = filt_leads{i}*1000;
            
            % Get difference for entire signal
            noise_diff = raw_signal - filt_signal;
            noise_var(i) = var(noise_diff);
        end
        
        % Store results for this patient
        rms_min_all(patient) = min(flat_rms);  % Minimum RMS from flat spots
        rms_var_all(patient) = sqrt(mean(noise_var));  % RMS of variances
    end
    
    % Optional: Display some statistics
    %fprintf('Statistics across %d patients:\n', num_patients);
    %fprintf('RMS min - Mean: %.4f, Std: %.4f\n', mean(rms_min_all), std(rms_min_all));
    %fprintf('RMS var - Mean: %.4f, Std: %.4f\n', mean(rms_var_all), std(rms_var_all));
end