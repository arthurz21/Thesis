function [jt_anterior, jt_lateral] = JT_ant_lat_finder(ecg, sample_time, fidpts)

        j_start = fidpts(3);
        t_peak = fidpts(4);
        
        %anterior
        lead_V1_jt_area = sample_time*trapz(ecg.V1(j_start:t_peak));
        lead_V2_jt_area = sample_time*trapz(ecg.V2(j_start:t_peak));
        lead_V3_jt_area = sample_time*trapz(ecg.V3(j_start:t_peak));
        lead_V4_jt_area = sample_time*trapz(ecg.V4(j_start:t_peak));

        jt_anterior = (lead_V1_jt_area + lead_V2_jt_area + lead_V3_jt_area + lead_V4_jt_area)/4;

        %lateral
        lead_I_jt_area = sample_time*trapz(ecg.I(j_start:t_peak));
        lead_aVL_jt_area = sample_time*trapz(ecg.avL(j_start:t_peak));
        lead_V5_jt_area = sample_time*trapz(ecg.V5(j_start:t_peak));
        lead_V6_jt_area = sample_time*trapz(ecg.V6(j_start:t_peak));

        jt_lateral = (lead_I_jt_area + lead_aVL_jt_area + lead_V5_jt_area + lead_V6_jt_area)/4;

