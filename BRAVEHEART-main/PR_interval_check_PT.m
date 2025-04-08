function [PR_interval] = PR_interval_check_PT(raw_ecg)
    qrs_positions = detectQRSComplexes(raw_ecg, 1000, false);
    t_positions = detectTWaves(raw_ecg, 1000, qrs_positions, false);
    pvc_labels = detectPVCs(raw_ecg, 1000, qrs_positions, false);
    afib_labels = detectAFIB(qrs_positions, 1000, pvc_labels, false);
    p_wave_present = pathologyCheck(pvc_labels, afib_labels, false);
    p_wave_candidates = detectNormalPWave(raw_ecg, qrs_positions, 1000, p_wave_present, false);
    diss_p_candidates = detectDissociatedPWave(raw_ecg, qrs_positions, 1000, p_wave_present, pvc_labels, false);
    final_p_waves = verifyPWaves(raw_ecg, qrs_positions, p_wave_candidates, 1000, false);
    PR_intervals = compute_PR_intervals(raw_ecg, final_p_waves, qrs_positions, 1000, false);
    
    
    if isempty(PR_intervals)
        PR_interval = 0; % If no valid PR intervals, return NaN
    else
        PR_interval = mean(PR_intervals);  % Calculate average PR interval
    end

end 


function [qrs_positions] = detectQRSComplexes(ecg_signal, fs, debug_plots)
    % QRS complex detection using Phasor Transform with improved R peak discrimination
    % Input:
    %   ecg_signal: Raw ECG signal
    %   fs: Sampling frequency (Hz)
    %   debug_plots: Boolean to show intermediate plots (default: true)
    % Output:
    %   qrs_positions: Positions of detected R peaks only
    
    if nargin < 3
        debug_plots = true;
    end
    
    % Make sure signal is a column vector
    if size(ecg_signal,2) > size(ecg_signal,1)
        ecg_signal = ecg_signal';
    end
    
    % Remove any DC offset
    ecg_signal = ecg_signal - mean(ecg_signal);
    
    % Normalize signal amplitude
    ecg_signal = ecg_signal / max(abs(ecg_signal));
    
    % Step 1: Preprocessing
    % High-pass filter to remove baseline wander (cutoff freq = 0.67 Hz)
    [b_high, a_high] = butter(4, 0.67/(fs/2), 'high');
    ecg_filtered = filtfilt(b_high, a_high, ecg_signal);
    
    % Bandpass filter focused more on QRS frequency range (10-25 Hz)
    [b_band, a_band] = butter(4, [10/(fs/2) 25/(fs/2)], 'bandpass');
    ecg_filtered = filtfilt(b_band, a_band, ecg_filtered);
    
    % Step 2: Apply Phasor Transform
    Rv = 0.001;
    phase_signal = zeros(size(ecg_filtered));
    phasor_complex = zeros(size(ecg_filtered));
    
    for n = 1:length(ecg_filtered)
        % Calculate complex phasor
        phasor_complex(n) = Rv + 1j * ecg_filtered(n);
        % Calculate phase angle
        phase_signal(n) = atan2(ecg_filtered(n), Rv);
    end
    
    % Calculate magnitude
    magnitude_signal = abs(phasor_complex);
    
    if debug_plots
        figure('Name', 'QRS Detection Steps');
        
        % Plot 1: Raw ECG
        subplot(4,1,1);
        plot((1:length(ecg_signal))/fs, ecg_signal);
        title('Raw ECG Signal');
        grid on;
        xlabel('Time (s)');
        ylabel('Amplitude');
        
        % Plot 2: Filtered ECG
        subplot(4,1,2);
        plot((1:length(ecg_filtered))/fs, ecg_filtered);
        title('Filtered ECG Signal');
        grid on;
        xlabel('Time (s)');
        ylabel('Amplitude');
        
        % Plot 3: Phasor Transform
        subplot(4,1,3);
        plot((1:length(phase_signal))/fs, phase_signal);
        title('Phase Signal after Phasor Transform');
        grid on;
        xlabel('Time (s)');
        ylabel('Phase (rad)');
        
        % Plot 4: Final Detection Results
        subplot(4,1,4);
        plot((1:length(ecg_signal))/fs, ecg_signal);
        hold on;
    end
    
    % Step 3: Find potential peaks with higher threshold
    window_length = round(0.2 * fs); % 200 ms window
    min_peak_distance = round(0.2 * fs); % Minimum 200ms between peaks
    
    % Use findpeaks with minimum peak distance and height requirements
    [peak_values, peak_locations] = findpeaks(ecg_filtered, ...
        'MinPeakDistance', min_peak_distance, ...
        'MinPeakHeight', 0.4 * max(ecg_filtered)); % Require peaks to be at least 40% of max amplitude
    
    % Further refine peaks based on relative amplitude
    median_peak = median(peak_values);
    qrs_positions = peak_locations(peak_values > 0.5 * max(peak_values));
    
    % Add detection results to the final plot
    if debug_plots
        if ~isempty(qrs_positions)
            plot(qrs_positions/fs, ecg_signal(qrs_positions), 'ro', 'MarkerSize', 10);
        end
        title('QRS Detection Results (R peaks only)');
        legend('ECG Signal', 'Detected R peaks');
        grid on;
        xlabel('Time (s)');
        ylabel('Amplitude');
    end
    
    % Calculate basic heart rate statistics if peaks were found
    if ~isempty(qrs_positions)
        RR_intervals = diff(qrs_positions) / fs;
        mean_HR = 60 / mean(RR_intervals);
        %fprintf('Mean heart rate: %.1f BPM\n', mean_HR);
        %fprintf('Number of R peaks detected: %d\n', length(qrs_positions));
    else
        %warning('No R peaks detected');
    end

end

function [t_positions] = detectTWaves(ecg_signal, fs, qrs_positions, debug_plots)
    % T wave detection using Phasor Transform
    % Input:
    %   ecg_signal: Raw ECG signal
    %   fs: Sampling frequency (Hz)
    %   qrs_positions: R peak positions from QRS detection
    %   debug_plots: Boolean to show plots (default: true)
    % Output:
    %   t_positions: Detected T wave peak positions
    
    if nargin < 4
        debug_plots = true;
    end
    
    % Make signal a column vector if needed
    if size(ecg_signal,2) > size(ecg_signal,1)
        ecg_signal = ecg_signal';
    end
    
    % Remove baseline wander
    [b_high, a_high] = butter(4, 0.67/(fs/2), 'high');
    ecg_filtered = filtfilt(b_high, a_high, ecg_signal);
    
    % Smooth signal with median filter
    window_size = round(0.04 * fs); % 40ms window
    ecg_filtered = medfilt1(ecg_filtered, window_size);
    
    % Initialize arrays
    t_positions = [];  % Changed to empty array instead of preallocating
    phase_signals = {};
    
    % Apply Phasor Transform for T wave detection
    Rv = 0.1; % Value specified in the paper for T waves
    
    % Process each QRS complex
    for i = 1:length(qrs_positions)
        % Calculate search window for current T wave
        if i < length(qrs_positions)
            RR = qrs_positions(i+1) - qrs_positions(i);
        else
            RR = mean(diff(qrs_positions)); % Use mean RR for last beat
        end
        
        % Define search window (paper's specifications)
        start_idx = round(qrs_positions(i) + 0.12 * RR);
        end_idx = min(length(ecg_signal), round(qrs_positions(i) + 0.57 * RR + 0.06 * fs));
        
        % Check if indices are valid
        if end_idx <= start_idx || start_idx < 1 || start_idx > length(ecg_signal) || end_idx > length(ecg_signal) || ~isfinite(start_idx) || ~isfinite(end_idx)
            continue;
        end
        
        % Extract window for T wave detection
        window = ecg_filtered(start_idx:end_idx);
        
        % Skip if window is empty
        if isempty(window)
            continue;
        end
        
        % Apply Phasor Transform to window
        phase_signal = zeros(size(window));
        for n = 1:length(window)
            phase_signal(n) = atan2(window(n), Rv);
        end
        
        if i <= length(phase_signals)
            phase_signals{i} = phase_signal;
        else
            phase_signals{end+1} = phase_signal;
        end
        
        % Find maximum in phase signal (T wave peak)
        [~, max_idx] = max(phase_signal);
        t_position = start_idx + max_idx - 1;
        
        % Add to results only if position is valid
        if t_position > 0 && t_position <= length(ecg_signal) && isfinite(t_position)
            t_positions = [t_positions; t_position];
        end
    end
    
    % Plotting
    if debug_plots
        figure('Name', 'T Wave Detection');
        
        % Plot 1: Full signal with detected waves
        subplot(2,1,1);
        plot((1:length(ecg_signal))/fs, ecg_signal);
        hold on;
        plot(qrs_positions/fs, ecg_signal(qrs_positions), 'ro', 'MarkerSize', 10);
        if ~isempty(t_positions)
            plot(t_positions/fs, ecg_signal(t_positions), 'go', 'MarkerSize', 10);
        end
        title('ECG with Detected QRS (red) and T waves (green)');
        xlabel('Time (s)');
        ylabel('Amplitude');
        grid on;
        legend('ECG', 'R peaks', 'T waves');
        
        % Plot 2: Example of T wave detection in a single beat
        if ~isempty(qrs_positions)
            subplot(2,1,2);
            % Take middle beat for demonstration
            beat_idx = round(length(qrs_positions)/2);
            if beat_idx >= 1 && beat_idx <= length(qrs_positions)
                if beat_idx < length(qrs_positions)
                    RR = qrs_positions(beat_idx+1) - qrs_positions(beat_idx);
                else
                    RR = mean(diff(qrs_positions));
                end
                start_idx = round(qrs_positions(beat_idx) + 0.12 * RR);
                end_idx = min(length(ecg_signal), round(qrs_positions(beat_idx) + 0.57 * RR + 0.06 * fs));
                
                if start_idx < end_idx && start_idx >= 1 && end_idx <= length(ecg_signal)
                    % Plot the window
                    t = (start_idx:end_idx)/fs;
                    plot(t, ecg_signal(start_idx:end_idx));
                    hold on;
                    
                    % Find this t_position
                    t_idx = [];
                    if ~isempty(t_positions)
                        t_idx = find(t_positions >= start_idx & t_positions <= end_idx, 1);
                        if ~isempty(t_idx)
                            plot(t_positions(t_idx)/fs, ecg_signal(t_positions(t_idx)), 'go', 'MarkerSize', 10);
                        end
                    end
                    
                    title('Example T Wave Detection Window');
                    xlabel('Time (s)');
                    ylabel('Amplitude');
                    grid on;
                end
            end
        end
    end
end

function [pvc_labels] = detectPVCs(ecg_signal, fs, qrs_positions, debug_plots)
    % PVC detection using Area Under Curve method
    % Input:
    %   ecg_signal: Raw ECG signal
    %   fs: Sampling frequency (Hz)
    %   qrs_positions: R peak positions from QRS detection
    %   debug_plots: Boolean to show plots (default: true)
    % Output:
    %   pvc_labels: Boolean array indicating PVC beats (true/false)
    
    if nargin < 4
        debug_plots = true;
    end
    
    % Initialize output
    pvc_labels = false(size(qrs_positions));
    
    % Remove baseline wander using high-pass filter
    [b_high, a_high] = butter(4, 0.67/(fs/2), 'high');
    ecg_filtered = filtfilt(b_high, a_high, ecg_signal);
    
    % Calculate window size in samples (±150ms)
    window_samples = round(0.15 * fs);
    
    % Calculate AUC for each QRS complex
    aucs = zeros(size(qrs_positions));
    for i = 1:length(qrs_positions)
        % Define window boundaries
        start_idx = max(1, qrs_positions(i) - window_samples);
        end_idx = min(length(ecg_filtered), qrs_positions(i) + window_samples);
        
        % Calculate AUC using trapezoidal integration
        segment = abs(ecg_filtered(start_idx:end_idx));
        aucs(i) = trapz(segment);
    end
    
    % Calculate adaptive threshold (1.3 * median AUC)
    auc_threshold = 1.3 * median(aucs);
    
    % Detect PVCs
    pvc_labels = aucs > auc_threshold;
    
    % Check if more than 75% of beats are labeled as PVCs
    % If so, consider it as bundle branch block and label all as normal
    if mean(pvc_labels) > 0.75
        pvc_labels(:) = false;
    end
    
    % Debug plotting
    if debug_plots
        figure('Name', 'PVC Detection');
        
        % Plot 1: Full signal with detected PVCs
        subplot(2,1,1);
        t = (1:length(ecg_signal))/fs;
        plot(t, ecg_signal);
        hold on;
        
        % Plot R peaks with color coding (red for PVC, green for normal)
        for i = 1:length(qrs_positions)
            if pvc_labels(i)
                plot(qrs_positions(i)/fs, ecg_signal(qrs_positions(i)), 'go', 'MarkerSize', 10);
            else
                plot(qrs_positions(i)/fs, ecg_signal(qrs_positions(i)), 'ro', 'MarkerSize', 10);
            end
        end
        
        title('ECG with PVC Detection');
        xlabel('Time (s)');
        ylabel('Amplitude');
        % Create legend with specific entries
        h_ecg = plot(nan, nan, 'b-', 'DisplayName', 'ECG');
        h_pvc = plot(nan, nan, 'go', 'MarkerSize', 10, 'DisplayName', 'PVC beats');
        h_normal = plot(nan, nan, 'ro', 'MarkerSize', 10, 'DisplayName', 'Normal beats');
        legend([h_ecg, h_pvc, h_normal]);
        grid on;
        
        % Plot 2: AUC values and threshold
        subplot(2,1,2);
        beat_nums = 1:length(aucs);
        plot(beat_nums, aucs, 'b.-');
        hold on;
        plot(beat_nums([1 end]), [auc_threshold auc_threshold], 'r--');
        title('Area Under Curve Analysis');
        xlabel('Beat Number');
        ylabel('AUC Value');
        legend('Beat AUC', 'PVC Threshold');
        grid on;
    end
end

function [afib_labels] = detectAFIB(qrs_positions, fs, pvc_labels, debug_plots)
    % AFIB detection using Shannon entropy of heart rate dynamics
    % Input:
    %   qrs_positions: R peak positions
    %   fs: Sampling frequency (Hz)
    %   pvc_labels: Boolean array indicating PVC beats
    %   debug_plots: Boolean to show plots (default: true)
    % Output:
    %   afib_labels: Boolean array indicating AFIB beats (true/false)
    
    if nargin < 4
        debug_plots = true;
    end
    
    % Initialize output
    afib_labels = false(size(qrs_positions));
    
    % Calculate RR intervals in seconds
    rr_intervals = diff(qrs_positions) / fs;
    
    % Calculate heart rate sequence
    hr = 60 ./ rr_intervals;
    
    % Create symbolic sequence using 3 symbols
    % Normalize heart rate to 0-1 range for symbolization
    hr_norm = (hr - min(hr)) / (max(hr) - min(hr));
    symbols = zeros(size(hr_norm));
    symbols(hr_norm <= 1/3) = 0;
    symbols((hr_norm > 1/3) & (hr_norm <= 2/3)) = 1;
    symbols(hr_norm > 2/3) = 2;
    
    % Calculate Shannon entropy in sliding window
    window_size = 59; % as specified in the paper
    sh_values = zeros(size(qrs_positions));
    
    for i = (window_size+1)/2:length(qrs_positions)-(window_size-1)/2
        % Get current window of symbols
        start_idx = i - (window_size-1)/2;
        end_idx = i + (window_size-1)/2;
        if end_idx > length(symbols)
            break;
        end
        window = symbols(start_idx:end_idx);
        
        % Calculate word distribution
        words = [];
        for j = 1:length(window)-2
            word = window(j:j+2);
            words = [words; word];
        end
        
        % Convert words to strings for unique counting
        word_strings = arrayfun(@(i) num2str(words(i,:)), 1:size(words,1), 'UniformOutput', false);
        
        % Calculate probabilities
        [unique_words, ~, ic] = unique(word_strings);
        word_counts = histcounts(ic, 1:length(unique_words)+1);
        probabilities = word_counts / sum(word_counts);
        
        % Calculate Shannon entropy
        sh = -sum(probabilities .* log2(probabilities + eps));
        sh_values(i) = sh;
    end
    
    % Mark AFIB where Shannon entropy exceeds threshold
    threshold = 0.737; % as specified in the paper
    afib_labels(sh_values > threshold) = true;
    
    % Check for PVC influence
    for i = (window_size+1)/2:length(qrs_positions)-(window_size-1)/2
        if afib_labels(i)
            % Count PVCs in current window
            start_idx = max(1, i - (window_size-1)/2);
            end_idx = min(length(pvc_labels), i + (window_size-1)/2);
            window_pvcs = sum(pvc_labels(start_idx:end_idx));
            
            % If more than 30 PVCs in window, consider it PVC influence not AFIB
            if window_pvcs > 30
                afib_labels(i) = false;
            end
        end
    end
    
    % Debug plotting
    if debug_plots
        figure('Name', 'AFIB Detection');
        
        % Plot 1: RR intervals
        subplot(3,1,1);
        t = cumsum([0; rr_intervals]);
        plot(t(1:end-1), rr_intervals, 'b.-');
        title('RR Intervals');
        xlabel('Time (s)');
        ylabel('RR Interval (s)');
        grid on;
        
        % Plot 2: Shannon entropy
        subplot(3,1,2);
        plot(1:length(sh_values), sh_values, 'b.-');
        hold on;
        plot([1 length(sh_values)], [threshold threshold], 'r--');
        title('Shannon Entropy');
        xlabel('Beat Number');
        ylabel('Entropy Value');
        legend('Shannon Entropy', 'AFIB Threshold');
        grid on;
        
        % Plot 3: AFIB detection result
        subplot(3,1,3);
        plot(t(1:end-1), double(afib_labels(1:end-1)), 'b.-');
        title('AFIB Detection Result');
        xlabel('Time (s)');
        ylabel('AFIB Present');
        ylim([-0.1 1.1]);
        grid on;
    end
end

function [p_wave_present] = pathologyCheck(pvc_labels, afib_labels, debug_plots)
    % Pathology check for determining if P wave should be present
    % Input:
    %   pvc_labels: Boolean array indicating PVC beats (true/false)
    %   afib_labels: Boolean array indicating AFIB beats (true/false)
    % Output:
    %   p_wave_present: Boolean array indicating if P wave should be present
    
    % Initialize output array
    p_wave_present = false(size(pvc_labels));
    
    % Check each beat
    for i = 1:length(pvc_labels)
        % According to the paper's flowchart:
        % 1. If beat is PVC, P wave is not present
        % 2. If beat is AFIB, P wave is not present
        % 3. Otherwise, P wave should be present
        
        if pvc_labels(i)
            p_wave_present(i) = false;  % No P wave during PVC
        elseif afib_labels(i)
            p_wave_present(i) = false;  % No P wave during AFIB
        else
            p_wave_present(i) = true;   % P wave should be present
        end
    end
    
    if debug_plots
        % Optional: Create visualization to show pathology check results
        figure('Name', 'Pathology Check Results');
        
        % Create time axis for x-axis
        t = 1:length(p_wave_present);
        
        % Plot pathology check results
        subplot(3,1,1);
        stem(t, pvc_labels, 'filled');
        title('PVC Labels');
        ylim([-0.1 1.1]);
        grid on;
        
        subplot(3,1,2);
        stem(t, afib_labels, 'filled');
        title('AFIB Labels');
        ylim([-0.1 1.1]);
        grid on;
        
        subplot(3,1,3);
        stem(t, p_wave_present, 'filled');
        title('P Wave Present');
        xlabel('Beat Number');
        ylim([-0.1 1.1]);
        grid on;
    end
end

function [p_wave_candidates] = detectNormalPWave(ecg, qrs_positions, fs, p_wave_present, debug_plots)
    % Normal P wave detection using phasor transform
    % Input:
    %   ecg: ECG signal
    %   qrs_positions: R peak positions
    %   fs: Sampling frequency (Hz)
    %   p_wave_present: Boolean array indicating where to look for P waves
    % Output:
    %   p_wave_candidates: Detected P wave positions

    % Initialize output array with same length as qrs_positions
    p_wave_candidates = [];
    
    % Check if there are enough QRS positions for meaningful analysis
    if length(qrs_positions) < 2
        return;  % Not enough QRS complexes to calculate RR intervals
    end
    
    % Calculate RR intervals
    rr_intervals = diff(qrs_positions);
    
    % Add first RR interval as a copy of the first one for the first beat
    rr_intervals = [rr_intervals(1); rr_intervals];
    
    % Parameters for phasor transform
    Rv = 0.05; % as specified in the paper
    
    % Process each beat
    for i = 2:length(qrs_positions)  % Start from second beat
        % Skip if P wave shouldn't be present
        if i > length(p_wave_present) || ~p_wave_present(i)
            continue;
        end
        
        % Area demarcation for P wave searching
        % R(i-1)+0.71×RR(i) to R(i)-0.07×RR(i)-60ms
        start_idx = round(qrs_positions(i-1) + 0.71 * rr_intervals(i));
        end_idx = round(qrs_positions(i) - 0.07 * rr_intervals(i) - 0.06 * fs);
        
        % Ensure indices are within signal bounds
        start_idx = max(1, min(length(ecg), start_idx));
        end_idx = max(1, min(length(ecg), end_idx));
        
        if end_idx <= start_idx
            continue;
        end
        
        % Extract search window
        search_window = ecg(start_idx:end_idx);
        
        % Apply phasor transform
        % y(n) = Rv + jx(n)
        y = Rv + 1j * search_window;
        
        % Calculate phase signal PT(n) = tan^(-1)(x(n)/Rv)
        phase = atan(search_window / Rv);
        
        % Find maximum in phase signal
        [~, max_idx] = max(phase);
        
        % Convert back to original signal index
        p_wave_candidates = [p_wave_candidates; start_idx + max_idx - 1];
    end
    
    % Remove zeros (no detection)
    p_wave_candidates = p_wave_candidates(p_wave_candidates > 0);
    
    if debug_plots && ~isempty(ecg)
        % Create time vector
        t = (0:length(ecg)-1)/fs;

        % Visualization
        figure('Name', 'Normal P Wave Detection');
        
        % Plot 1: Original ECG with search windows
        subplot(2,1,1);
        plot(t, ecg, 'b', 'LineWidth', 1);
        hold on;
        
        % Plot QRS positions
        plot(t(qrs_positions), ecg(qrs_positions), 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        
        % Plot P wave candidates
        if ~isempty(p_wave_candidates)
            plot(t(p_wave_candidates), ecg(p_wave_candidates), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
        end
        
        % Highlight search windows
        for i = 2:length(qrs_positions)
            if i <= length(p_wave_present) && p_wave_present(i)
                start_idx = round(qrs_positions(i-1) + 0.71 * rr_intervals(i));
                end_idx = round(qrs_positions(i) - 0.07 * rr_intervals(i) - 0.06 * fs);
                
                % Ensure indices are within bounds
                start_idx = max(1, min(length(ecg), start_idx));
                end_idx = max(1, min(length(ecg), end_idx));
                
                if start_idx < end_idx
                    % Create patch coordinates
                    y_min = min(ecg);
                    y_max = max(ecg);
                    patch([t(start_idx) t(end_idx) t(end_idx) t(start_idx)], ...
                          [y_min y_min y_max y_max], ...
                          'green', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
            end
        end
        
        title('ECG Signal with P Wave Search Windows');
        xlabel('Time (s)');
        ylabel('Amplitude');
        legend('ECG', 'QRS', 'P waves', 'Search windows', 'Location', 'best');
        grid on;
        xlim([t(1) t(end)]);
        
        % Plot 2: PR intervals if any P waves were detected
        subplot(2,1,2);
        % Calculate PR intervals only for beats where P waves were detected
        pr_intervals = [];
        beat_numbers = [];
        
        for i = 1:length(p_wave_candidates)
            % Find corresponding QRS peak for this P wave
            next_qrs = qrs_positions(qrs_positions > p_wave_candidates(i));
            if ~isempty(next_qrs)
                pr_interval = (next_qrs(1) - p_wave_candidates(i)) / fs;
                pr_intervals = [pr_intervals; pr_interval];
                % Find which beat number this is
                beat_num = find(qrs_positions == next_qrs(1));
                beat_numbers = [beat_numbers; beat_num];
            end
        end
        
        % Plot PR intervals
        if ~isempty(pr_intervals) && ~isempty(beat_numbers)
            plot(beat_numbers, pr_intervals, 'bo-', 'LineWidth', 1, 'MarkerSize', 6);
            title('PR Intervals Over Time');
            xlabel('Beat Number');
            ylabel('PR Interval (seconds)');
            grid on;
            ylim([0 max(pr_intervals)*1.2]);
        else
            text(0.5, 0.5, 'No PR intervals detected', 'HorizontalAlignment', 'center');
            axis off;
        end
    end 
end

function [diss_p_candidates] = detectDissociatedPWave(ecg, qrs_positions, fs, p_wave_present, pvc_labels, debug_plots)
    % Initialize output
    diss_p_candidates = [];
    
    % Check if there are enough QRS positions for meaningful analysis
    if length(qrs_positions) < 2
        return;  % Not enough QRS complexes to calculate RR intervals
    end
    
    % Calculate RR intervals and handle first interval
    rr_intervals = diff(qrs_positions) / fs;  % Convert to seconds
    
    % Add first RR interval as a copy of the first one for the first beat
    if ~isempty(rr_intervals)
        rr_intervals = [rr_intervals(1); rr_intervals];  % First interval same as second
    else
        return; % Exit function if rr_intervals is empty
    end
    
    % Parameters for PT
    Rv = 0.05;
    
    % Process each RR interval starting from second beat
    for i = 2:(length(qrs_positions)-1)  % Stop one before last to avoid indexing error
        % Check three criteria for dissociated P waves
        if rr_intervals(i) > 1.6 * rr_intervals(i-1) && ...
           rr_intervals(i) > 1.6 && ...
           ~pvc_labels(i)
            
            % In search window: T(i-1)+200ms to P(i)-400ms
            start_idx = qrs_positions(i-1) + round(0.3 * fs);  % Approximate T wave location
            end_idx = qrs_positions(i) - round(0.4 * fs);
            
            if end_idx > start_idx
                % Extract window and apply PT
                window = ecg(start_idx:end_idx);
                phase = atan(window / Rv);
                [~, max_idx] = max(phase);
                diss_p_candidates = [diss_p_candidates; start_idx + max_idx - 1];
            end
        end
    end
    if debug_plots
        % Add visualization
        figure;
        t = (0:length(ecg)-1)/fs;
        plot(t, ecg, 'b');
        hold on;
        plot(t(qrs_positions), ecg(qrs_positions), 'rv', 'MarkerFaceColor', 'r');
        h = plot(nan, nan, 'mo', 'MarkerFaceColor', 'm');  % Create invisible plot for legend
        if ~isempty(diss_p_candidates)
            plot(t(diss_p_candidates), ecg(diss_p_candidates), 'mo', 'MarkerFaceColor', 'm');
        end
        title('ECG with Dissociated P Waves');
        xlabel('Time (s)');
        ylabel('Amplitude');
        legend('ECG', 'QRS', 'Dissociated P');
        grid on;
    end
end

function [final_p_waves] = verifyPWaves(ecg, qrs_positions, p_candidates, fs, debug_plots)
    % Initialize output
    final_p_waves = [];
    
    % Process each P wave candidate
    for i = 1:length(p_candidates)
        % Find corresponding QRS peak
        next_qrs = qrs_positions(qrs_positions > p_candidates(i));
        if isempty(next_qrs)
            continue;
        end
        qrs_peak = next_qrs(1);
        
        % Check amplitude criterion
        p_amp = abs(ecg(p_candidates(i)));
        qrs_amp = abs(ecg(qrs_peak));
        
        if p_amp > 0.05 * qrs_amp  % First criterion from paper
            % For now, we'll skip T wave check since we don't have T positions
            % In practice, would add: if P(i) > T(i-1)
            final_p_waves = [final_p_waves; p_candidates(i)];
        end
    end
    if debug_plots
        % Visualization
        figure;
        t = (0:length(ecg)-1)/fs;
        plot(t, ecg, 'b');
        hold on;
        plot(t(qrs_positions), ecg(qrs_positions), 'rv', 'MarkerFaceColor', 'r');
        plot(t(p_candidates), ecg(p_candidates), 'go', 'MarkerFaceColor', 'none', 'DisplayName', 'P candidates');
        plot(t(final_p_waves), ecg(final_p_waves), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Verified P');
        title('ECG with Verified P Waves');
        xlabel('Time (s)');
        ylabel('Amplitude');
        legend('ECG', 'QRS', 'P candidates', 'Verified P');
        grid on;
    end
end

function PR_intervals = compute_PR_intervals(ecg_signal, final_p_waves, qrs_positions, fs, debug_plots)
    % Initialize arrays for onsets
    p_wave_onsets = zeros(size(final_p_waves));
    qrs_onsets = zeros(size(qrs_positions));
    PR_intervals = [];
    
    % Parameters
    window_before_p = round(0.05 * fs);  % 50ms window before P peak
    window_before_r = round(0.08 * fs);  % 80ms window before R peak (increased to catch Q waves)
    
    % 1. Detect P wave onsets (keeping original code)
    for i = 1:length(final_p_waves)
        p_peak = final_p_waves(i);
        start_idx = max(1, p_peak - window_before_p);
        
        if start_idx < p_peak
            window = ecg_signal(start_idx:p_peak);
            d_window = diff(window);
            threshold = 0.1 * max(abs(d_window));
            onset_idx = find(abs(d_window) > threshold, 1, 'first');
            
            if ~isempty(onset_idx)
                p_wave_onsets(i) = start_idx + onset_idx - 1;
            else
                p_wave_onsets(i) = start_idx;
            end
        end
    end
    
    % 2. Detect QRS onsets with improved Q wave detection
    for i = 1:length(qrs_positions)
        r_peak = qrs_positions(i);
        start_idx = max(1, r_peak - window_before_r);
        
        if start_idx < r_peak
            % Extract window
            window = ecg_signal(start_idx:r_peak);
            
            % Calculate first derivative
            d_window = diff(window);
            
            % Look for significant negative deflection first (Q wave)
            neg_threshold = -0.15 * max(abs(d_window));  % Threshold for negative deflection
            q_candidate_idx = find(d_window < neg_threshold, 1, 'first');
            
            if ~isempty(q_candidate_idx)
                % Q wave found - this is our onset
                qrs_onsets(i) = start_idx + q_candidate_idx - 1;
            else
                % No clear Q wave, look for any significant deflection
                threshold = 0.2 * max(abs(d_window));
                onset_idx = find(abs(d_window) > threshold, 1, 'first');
                
                if ~isempty(onset_idx)
                    qrs_onsets(i) = start_idx + onset_idx - 1;
                else
                    % If no clear onset found, look for the point where slope starts changing
                    slope_changes = diff(sign(d_window));
                    slope_change_idx = find(slope_changes ~= 0, 1, 'first');
                    
                    if ~isempty(slope_change_idx)
                        qrs_onsets(i) = start_idx + slope_change_idx - 1;
                    else
                        qrs_onsets(i) = start_idx;
                    end
                end
            end
        end
    end
    
    % 3. Calculate PR intervals (keeping original code)
    for i = 1:length(p_wave_onsets)
        next_qrs_onsets = qrs_onsets(qrs_onsets > p_wave_onsets(i));
        
        if ~isempty(next_qrs_onsets)
            pr_interval = (next_qrs_onsets(1) - p_wave_onsets(i)) / fs * 1000;
            
            if pr_interval > 0
                PR_intervals = [PR_intervals; pr_interval];
            end
        end
    end
    
    % Debug plotting (keeping original code)
    if debug_plots
        figure('Name', 'PR Interval Detection');
        t = (1:length(ecg_signal))/fs;
        plot(t, ecg_signal, 'b');
        hold on;
        
        plot(final_p_waves/fs, ecg_signal(final_p_waves), 'go', 'DisplayName', 'P Peaks');
        plot(p_wave_onsets/fs, ecg_signal(p_wave_onsets), 'g*', 'DisplayName', 'P Onsets');
        plot(qrs_positions/fs, ecg_signal(qrs_positions), 'ro', 'DisplayName', 'R Peaks');
        plot(qrs_onsets/fs, ecg_signal(qrs_onsets), 'r*', 'DisplayName', 'QRS Onsets');
        
        title('ECG with PR Interval Detection Points');
        xlabel('Time (s)');
        ylabel('Amplitude');
        legend('Location', 'best');
        grid on;
        
        for i = 1:min(3, length(PR_intervals))
            x_start = p_wave_onsets(i)/fs;
            x_end = qrs_onsets(find(qrs_onsets > p_wave_onsets(i), 1))/fs;
            y_pos = min(ecg_signal) - 0.1*range(ecg_signal);
            
            plot([x_start x_end], [y_pos y_pos], 'k-');
            text((x_start + x_end)/2, y_pos, sprintf('%.0f ms', PR_intervals(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        end
    end
end