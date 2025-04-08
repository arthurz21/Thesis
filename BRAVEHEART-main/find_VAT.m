function [vat_II, vat_III, vat_V2, vat_V4] = find_VAT(ecg, QRST_pts)
    q_start = QRST_pts(1);
    s_end = QRST_pts(3);

    max_ecg_II  = max(ecg.II(q_start:s_end));
    max_ecg_III = max(ecg.III(q_start:s_end));
    max_ecg_V2 = max(ecg.V2(q_start:s_end));
    max_ecg_V4 = max(ecg.V4(q_start:s_end));
    
    R_peak_loc_II = find(ecg.II == max_ecg_II);
    R_peak_loc_III = find(ecg.III == max_ecg_III);
    R_peak_loc_V2 = find(ecg.V2 == max_ecg_V2);
    R_peak_loc_V4 = find(ecg.V4 == max_ecg_V4);



    vat_II = (R_peak_loc_II - q_start)*(1000/(ecg.hz));
    vat_III = (R_peak_loc_III - q_start)*(1000/(ecg.hz));
    vat_V2 = (R_peak_loc_V2 - q_start)*(1000/(ecg.hz));
    vat_V4 = (R_peak_loc_V4 - q_start)*(1000/(ecg.hz));