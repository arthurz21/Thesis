function [latConcaveAmp, antConcaveAmp] = computeConcaveAmp(ecg_leads, fidpts)
    % Inputs:
    % - ecg_leads: A struct with ECG signals for leads (e.g., I, aVL, V5, V6)
    % - fidpts: A struct with indices for Q start, S end, T start, T end, etc.
    % Output:
    % - latConcaveAmp: Sum of concave amplitudes in lateral leads

    % Define lateral leads
    lateral_leads = {'I', 'avL', 'V5', 'V6'};
    anterior_leads = {'V1', 'V2', 'V3', 'V4'};

    % Initialize sum
    latConcaveAmp = 0;
    antConcaveAmp = 0;

    % Loop through each lateral lead
    for lead = lateral_leads
        % Extract the signal for the current lead
        signal = ecg_leads.(lead{1});
        
        % Extract region of interest (e.g., Q start to S end, T start to T end)
        qrst_region = signal(fidpts(1):fidpts(5));
        
        % Compute second derivative
        second_derivative = diff(qrst_region, 2);
        
        % Find positive second derivative (concave regions)
        concave_indices = find(second_derivative > 0);

        % Sum the minimum amplitudes in the concave regions
        if ~isempty(concave_indices)
            concave_amplitudes = qrst_region(concave_indices+2);
            latConcaveAmp = latConcaveAmp + sum(concave_amplitudes);
        end
    end

    for lead = anterior_leads
        % Extract the signal for the current lead
        signal = ecg_leads.(lead{1});
        
        % Extract region of interest (e.g., Q start to S end, T start to T end)
        qrst_region = signal(fidpts(1):fidpts(5));
        
        % Compute second derivative
        second_derivative = diff(qrst_region, 2);
        
        % Find positive second derivative (concave regions)
        concave_indices = find(second_derivative > 0);
       

        % Sum the minimum amplitudes in the concave regions
        if ~isempty(concave_indices)
            concave_amplitudes = qrst_region(concave_indices + 2 );
            antConcaveAmp = latConcaveAmp + sum(concave_amplitudes); % Take absolute amplitudes
        end
    end

end
