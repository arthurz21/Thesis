function [rms_min, rms_var] = calculateECGMetrics(ecg, Q_onset, T_end)
    % Calculate RMS min and variance using 8 independent leads
    % Only considers active ECG segment between Q onset and T end
    %
    % Input:
    %   ecg - Median beat object with leads accessed as ecg.I, ecg.II, etc.
    %   Q_onset - Sample index where Q wave begins
    %   T_end - Sample index where T wave ends
    %
    % Output:
    %   rms_min - Minimum RMS value in the active segment
    %   rms_var - Variance of RMS values in the active segment
    
    % Get number of time points
    num_samples = length(ecg.I);  % assuming all leads have same length
    
    % Get all 8 independent leads as a matrix (samples × 8)
    independent_leads = [
        ecg.I, ecg.II, ...    % Limb leads
        ecg.V1, ecg.V2, ecg.V3, ecg.V4, ecg.V5, ecg.V6  % Precordial leads
    ]*1000;
    
    % Extract active segment
    active_leads = independent_leads(Q_onset:T_end, :);
    
    % Calculate RMS for each time point in active segment
    num_active_samples = size(active_leads, 1);
    rms_values = zeros(num_active_samples, 1);
    
    for t = 1:num_active_samples
        % Get values from all leads at time t
        values_at_t = active_leads(t, :);
        
        % Calculate RMS according to the formula
        rms_values(t) = sqrt(sum(values_at_t.^2) / 8);
    end
    
    % Calculate metrics
    rms_min = min(rms_values);
    rms_var = var(rms_values);
    
    % Optional: Display results
    %fprintf('RMS minimum: %.4f\n', rms_min);
    %fprintf('RMS variance: %.4f\n', rms_var);
    %fprintf('RMS values: %.4f\n', rms_values);
end