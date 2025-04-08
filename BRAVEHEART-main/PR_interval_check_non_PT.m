function [PR_interval] = PR_interval_check_non_PT(ecg, QRST_pts)
    %use first lead data to check for PR, as suggested by that
    %paper using phase transform
    QRS_start = QRST_pts(1);
    QRS_end = QRST_pts(3);
    R_peak = max(ecg.I(QRS_start:QRS_end));

    %now we check after T end, what is the maximum amplitude
    %compare 5% to the QRS peak to see if we can classify this as P peak
    T_end = QRST_pts(5);
    [P_peak, P_peak_index] = max(ecg.I(T_end:end));
    abs_P_peak_loc = T_end + P_peak_index - 1;

    if P_peak >= 0.05*R_peak
        PR_interval = QRS_start + length(ecg.I) - abs_P_peak_loc;
    else
        PR_interval = 0;
    end

end