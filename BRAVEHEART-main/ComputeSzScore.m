function [MIsz] = ComputeSzScore(ecg, QRST_pts, fidpts, isMale, age)
    %fprintf("checkpoint1\n")
    terminal_deflection = analyzeV1WavePattern(ecg, QRST_pts);
    
    %isMale = 0;
    RBBB_status = 0;
    LAFB_status = 0;
    LVH_status = 0;
    No_confounders = 0;
    LBBB_status = 0;
    RAO_status = 0;

    if terminal_deflection == 0
        LBBB_status = check_LBBB(ecg, QRST_pts, isMale);
        
        if LBBB_status == 0
            LAFB_status = check_LAFB(ecg, QRST_pts, isMale, fidpts);
        end

        if LAFB_status == 0
            LVH_status = check_LVH(ecg, isMale, fidpts);
        end

        if LVH_status == 0
            No_confounders = 1;
        else
            No_confounders = 0;
        end

    elseif terminal_deflection == 1
       QRS_duration = QRST_pts(3) - QRST_pts(1);
       
       if QRS_duration >= 120
           %Lead 1
           L1_signal = ecg.I;
           
           [L1_r_wave, L1_s_wave, ~, ~, ~] = rs_values(L1_signal,fidpts);
        
           %Lead avF
           avF_signal = ecg.avF;
           [avF_r_wave, avF_s_wave, ~, ~, ~] = rs_values(avF_signal,fidpts);

           qrs_frontal_axis = ecg_axis(L1_r_wave, L1_s_wave, avF_r_wave, avF_s_wave);
           
           if qrs_frontal_axis <= 45 && qrs_frontal_axis > 180
               RBBB_status = 1;
               LAFB_status = 1;
           else
               RBBB_status = 1;
           end

       else
            LAFB_status = check_LAFB(ecg, QRST_pts, isMale, fidpts);

            if LAFB_status == 0
                LVH_status = check_LVH(ecg, isMale, fidpts);
            end

            if LVH_status == 0
                No_confounders = 1;
            else
                No_confounders = 0;
            end
       end
     end
    
    RAO_status = check_RAO(ecg, QRST_pts);
    
    %put it in another function since this one is getting long
    
    MIsz = compute_score(ecg, QRST_pts, fidpts, isMale, RBBB_status, LAFB_status, LVH_status,...
        No_confounders, LBBB_status, RAO_status, age);
end

function [score] = compute_score(ecg, QRST_pts, fidpts, isMale, RBBB_status, LAFB_status, LVH_status,...
        No_confounders, LBBB_status, RAO_status, age)
    
    %remember RAO's effect for the criterias
    d_adj = (1 - 0.1 * (~isMale));
    ageAdjustment = (age < 20) + (age >= 20 & age <= 54) .* (1 + (age - 20) * 0.01) + (age > 54) .* (1 + 0.34 - (age - 54) * 0.01);
    a_adj = ageAdjustment .* (isMale + ~isMale * 0.9);

    %Lead I
    L1_points = 0;
    
    %check do we actually have Q wave first
    [Q_exist] = Q_wave_detection(ecg.I, QRST_pts);
    
    if Q_exist ~= 0
        
        [Q_duration, Q_peak] = extractQDuration(ecg.I, QRST_pts, 5);
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.I, QRST_pts, Q_duration, RBBB_status);
    
    if R_exist && Q_peak
        R_Q = R_peak/Q_peak;
    else
        R_Q = 0;
    end
    
    if RBBB_status == 1 && LAFB_status == 0
       if Q_duration >= 30*d_adj
            L1_points = L1_points + 1;
       end

       if R_Q <=1
           L1_points = L1_points + 1;
       elseif R_peak <= 0.2*a_adj   
           %we only check this if above statement is not correct, since we
           %only do one point calculations per box
           L1_points = L1_points + 1;
       end
 
    elseif LAFB_status == 1 && RBBB_status == 0
       if Q_duration >= 30*d_adj    
            L1_points = L1_points + 1;
       end

       if R_Q <=1
           L1_points = L1_points + 1;
       elseif R_peak <= 0.2*a_adj
           %we only check this if above statement is not correct, since we
           %only do one point calculations per box
           L1_points = L1_points + 1;
       end

    elseif LAFB_status == 1 && RBBB_status == 1
       if Q_duration >= 30*d_adj
            L1_points = L1_points + 1;
       end

       if R_Q <=1
           L1_points = L1_points + 1;
       elseif R_peak <= 0.2*a_adj
           %we only check this if above statement is not correct, since we
           %only do one point calculations per box
           L1_points = L1_points + 1;
       end

    elseif LVH_status == 1
       if Q_duration >= 30*d_adj
            L1_points = L1_points + 1;
       end

       if R_Q <=1
           L1_points = L1_points + 1;
       elseif R_peak <= 0.2*a_adj   
           %we only check this if above statement is not correct, since we
           %only do one point calculations per box
           L1_points = L1_points + 1;
       end

    elseif No_confounders == 1
       if Q_duration >= 30*d_adj
            L1_points = L1_points + 1;
       end

       if R_Q <=1
           L1_points = L1_points + 1;
       elseif R_peak <= 0.2*a_adj
           %we only check this if above statement is not correct, since we
           %only do one point calculations per box
           L1_points = L1_points + 1;
       end


    elseif LBBB_status == 1
        [Q_peak, ~, R_peak, ~, ~, S_peak, ~, ~] = getLBBBMeasurements(ecg.I, QRST_pts, 5);
        if ~Q_peak 
            L1_points = L1_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_peak ~= 0 && R_peak/S_peak <= 1)
            L1_points = L1_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 1.5) || (S_peak ~= 0 && R_peak/S_peak <= 1.5)
            L1_points = L1_points + 1;
        end
    else
        L1_points = 0;
    end
    fprintf('checkpoint here Lead I\n');

    %Lead II
    L2_points = 0;
    [Q_exist] = Q_wave_detection(ecg.II, QRST_pts);
    
    if Q_exist ~= 0
        [Q_duration, Q_peak] = extractQDuration(ecg.II, QRST_pts, 5);
        
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    if RBBB_status == 1 && LAFB_status == 0
       if Q_duration >= 40*d_adj
           L2_points = L2_points + 2;
       elseif Q_duration >= 30*d_adj
           L2_points = L2_points + 1;
       end
     
    elseif LAFB_status == 1 && RBBB_status == 0
       if Q_duration >= 40*d_adj    
           L2_points = L2_points + 2;
       elseif Q_duration >= 30*d_adj
           L2_points = L2_points + 1;
       end

    elseif LAFB_status == 1 && RBBB_status == 1
       if Q_duration >= 40*d_adj    
           L2_points = L2_points + 2;
       elseif Q_duration >= 30*d_adj    
           L2_points = L2_points + 1;
       end
    elseif LVH_status == 1
       if Q_duration >= 40*d_adj    
           L2_points = L2_points + 2;
       elseif Q_duration >= 30*d_adj
           L2_points = L2_points + 1;
       end    
    elseif No_confounders == 1
       if Q_duration >= 40*d_adj
           L2_points = L2_points + 2;
       elseif Q_duration >= 30*d_adj
           L2_points = L2_points + 1;
       end
     
    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.II, QRST_pts, 5);
        if Q_duration>= 40*d_adj
            L2_points = L2_points + 2;
            
        elseif Q_duration >= 30*d_adj
            L2_points = L2_points + 1;
        end
        if (Q_peak ~= 0 && abs(Q_peak) > 0 && R_peak/abs(Q_peak) <= 0.5) || ...
           (S_peak ~= 0 && abs(S_peak) > 0 && R_peak/abs(S_peak) <= 0.5)
            L2_points = L2_points + 1;
        end
    else
        L2_points = 0;
    end
    
    fprintf('checkpoint here Lead II\n');
    %Lead avL
    avL_points = 0;
    [Q_exist] = Q_wave_detection(ecg.avL, QRST_pts);
    
    if Q_exist ~= 0
        
        [Q_duration, Q_peak] = extractQDuration(ecg.avL, QRST_pts, 5);
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.avL, QRST_pts, Q_duration, RBBB_status);
    
    if R_exist && Q_peak
        R_Q = R_peak/Q_peak;

    else
        R_Q = 0;
    end
    if RBBB_status == 1 && LAFB_status == 0
        if Q_duration >= 30*d_adj   
           avL_points = avL_points + 1;    
        end

        if R_Q <= 1
            avL_points = avL_points + 1;   
        end

    elseif LAFB_status == 1 && RBBB_status == 0
        if Q_duration >= 40*d_adj
           avL_points = avL_points + 1;    
        end

        if R_Q <= 1
            avL_points = avL_points + 1;   
        end

    elseif LAFB_status == 1 && RBBB_status == 1
        if Q_duration >= 40*d_adj
           avL_points = avL_points + 1;    
        end

        if R_Q <= 1
            avL_points = avL_points + 1;   
        end

    elseif LVH_status == 1
        if Q_duration >= 40*d_adj
           avL_points = avL_points + 1;    
        end

        if R_Q <= 1
            avL_points = avL_points + 1;   
        end  

    elseif No_confounders == 1
        if Q_duration >= 30*d_adj
           avL_points = avL_points + 1;    
        end

        if R_Q <= 1
            avL_points = avL_points + 1;   
        end

    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.avL, QRST_pts, 5)
        if Q_duration>= 50*d_adj    
            avL_points = avL_points + 2;

        elseif Q_duration >= 40*d_adj
            avL_points = avL_points + 1;
        end
        if (Q_peak ~= 0 && abs(Q_peak) > 0 && R_peak/abs(Q_peak) <= 0.5) || ...
           (S_peak ~= 0 && abs(S_peak) > 0 && R_peak/abs(S_peak) <= 0.5)
            avL_points = avL_points + 2;

        elseif (Q_peak ~= 0  && R_peak/abs(Q_peak) <= 1) || ...
           (S_peak ~= 0  && R_peak/abs(S_peak) <= 1)
            avL_points = avL_points + 1;
        end
    else
        avL_points = 0;
    end
    fprintf('checkpoint here Lead avL\n');

    %Lead avF
    avF_points = 0;
    [Q_exist] = Q_wave_detection(ecg.avF, QRST_pts);
    
    if Q_exist ~= 0
        
        [Q_duration, Q_peak] = extractQDuration(ecg.avF, QRST_pts, 5);
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.avF, QRST_pts, Q_duration, RBBB_status);
    
    if R_exist && Q_peak
        R_Q = R_peak/Q_peak;

    else
        R_Q = 0;
    end

    if RBBB_status == 1 && LAFB_status == 0
        if Q_duration >= 50*d_adj
           avF_points = avF_points + 3;    
        elseif Q_duration >= 40*d_adj
            avF_points = avF_points + 2;
        elseif Q_duration >= 30*d_adj
            avF_points = avF_points + 1;
        end

        if R_Q <= 1
            avF_points = avF_points + 2;   
        elseif R_Q <= 2
            avF_points = avF_points + 1;
        end 

    elseif LAFB_status == 1 && RBBB_status == 0
        if Q_duration >= 50*d_adj
           avF_points = avF_points + 3;    
        elseif Q_duration >= 40*d_adj
            avF_points = avF_points + 2;
        elseif Q_duration >= 30*d_adj
            avF_points = avF_points + 1;
        end

        if R_Q <= 1
            avF_points = avF_points + 2;   
        elseif R_Q <= 2
            avF_points = avF_points + 1;
        end 

    elseif LAFB_status == 1 && RBBB_status == 1
        if Q_duration >= 50*d_adj
           avF_points = avF_points + 3;    
        elseif Q_duration >= 40*d_adj
            avF_points = avF_points + 2;
        elseif Q_duration >= 30*d_adj
            avF_points = avF_points + 1;
        end

        if R_Q <= 1
            avF_points = avF_points + 2;   
        elseif R_Q <= 2
            avF_points = avF_points + 1;
        end 

    elseif LVH_status == 1
        if Q_duration >= 60*d_adj
           avF_points = avF_points + 3;    
        elseif Q_duration >= 50*d_adj
            avF_points = avF_points + 2;
        elseif Q_duration >= 40*d_adj
            avF_points = avF_points + 1;
        end

        if R_Q <= 1
            avF_points = avF_points + 2;   
        elseif R_Q <= 2
            avF_points = avF_points + 1;
        end 

    elseif No_confounders == 1
        if Q_duration >= 50*d_adj
           avF_points = avF_points + 3;    
        elseif Q_duration >= 40*d_adj
            avF_points = avF_points + 2;
        elseif Q_duration >= 30*d_adj
            avF_points = avF_points + 1;
        end

        if R_Q <= 1
            avF_points = avF_points + 2;   
        elseif R_Q <= 2
            avF_points = avF_points + 1;
        end 

    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.avF, QRST_pts, 5);
        if Q_duration>= 50*d_adj
            avF_points = avF_points + 2;

        elseif Q_duration >= 40*d_adj
            avF_points = avF_points + 1;
        end

        if (Q_peak ~= 0 && abs(Q_peak) > 0 && R_peak/abs(Q_peak) <= 0.5) || ...
           (S_peak ~= 0 && abs(S_peak) > 0 && R_peak/abs(S_peak) <= 0.5)
            avF_points = avF_points + 1;
        end   
    else
        avF_points = 0;
    end
    fprintf('checkpoint here Lead avF\n');

    %Lead V1 Ant
    V1_ant_points = 0;
    [Q_exist] = Q_wave_detection(ecg.V1, QRST_pts);
    
    if Q_exist ~= 0
        
        [Q_duration, Q_peak] = extractQDuration(ecg.V1, QRST_pts, 5);
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.V1, QRST_pts, Q_duration, RBBB_status);
    
    %check for any Q
    %already been done with Q_exist
    %check for init R
    [R_duration, initR_duration, initR_peak_amp, initR_peak_loc] = detectRWave(ecg.V1, QRST_pts, Q_duration, 5);
    [hasNtchInit40] = detectNtchInit40(ecg.V1, QRST_pts, 5);

    if RBBB_status == 1 && LAFB_status == 0
        if Q_duration >= 50*d_adj   
            V1_ant_points = V1_ant_points + 2;
        end
        if Q_exist == 1
            V1_ant_points = V1_ant_points + 1;
        elseif initR_duration <= 20*d_adj
            V1_ant_points = V1_ant_points + 1;
        end
        
    elseif LAFB_status == 1 && RBBB_status == 0
        if Q_exist == 1 && R_exist == 1
            V1_ant_points = V1_ant_points + 1;
        end

    elseif LAFB_status == 1 && RBBB_status == 1
        if Q_duration > 50*d_adj
            V1_ant_points = V1_ant_points + 2;
        end
        if Q_exist == 1
            V1_ant_points = V1_ant_points + 1;
        end

    elseif LVH_status == 1
        if Q_exist == 1 && R_exist == 1
            V1_ant_points = V1_ant_points + 1;
        elseif Q_exist == 1
            V1_ant_points = V1_ant_points + 1;
        elseif hasNtchInit40 == 1
            V1_ant_points = V1_ant_points + 1;

        end

    elseif No_confounders == 1
        if Q_exist == 1
            V1_ant_points = V1_ant_points + 1;
        end
    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.V1, QRST_pts, 5);
        if hasNchInit40 == 1
            V1_ant_points = V1_ant_points + 1;
        end
        if R_duration>= 30*d_adj || R_peak >= 0.3*a_adj
            V1_ant_points = V1_ant_points + 2;

        elseif R_duration>= 20*d_adj || R_peak >= 0.2*a_adj
            V1_ant_points = V1_ant_points + 1;
   
        end

    else
        V1_ant_points = 0;
    end
    
    fprintf('checkpoint here Lead V1 ant\n');
    %Lead V1 post
    %will reuse 
    V1_post_points = 0;
    [S_amplitude] = detectSWave(ecg.V1, QRST_pts, initR_peak_loc);
    
    if RAO_status == 1
        a = 0;
    elseif RBBB_status == 1 && LAFB_status == 0
        if initR_duration >= 60*d_adj || initR_peak_amp >= 1.5*a_adj
            V1_post_points = V1_post_points+2;

        elseif initR_peak_amp >= 1.0*d_adj || initR_duration >= 50*a_adj
            V1_post_points = V1_post_points+2;
        end
        
    elseif LAFB_status == 1 && RBBB_status == 0
        if (abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude >= 1)
            V1_post_points = V1_post_points+1;
        end
        if R_duration >= 50*d_adj || R_peak >= 1*a_adj
            V1_post_points = V1_post_points+2;
        elseif R_duration >= 40*d_adj || R_peak >= 0.7*a_adj
            V1_post_points = V1_post_points+1;
        end
        if Q_peak <= 0.2*a_adj && S_amplitude <= 0.2*a_adj
            V1_post_points = V1_post_points+1;
        end
     
    elseif LAFB_status == 1 && RBBB_status == 1
        if initR_duration >= 60*d_adj || initR_peak_amp >= 1.5*a_adj
            V1_post_points = V1_post_points+2;

        elseif initR_duration >= 50*d_adj || initR_peak_amp >= 1.0*a_adj
            V1_post_points = V1_post_points+1;
        end

    elseif LVH_status == 1
        if (abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude >= 1)
            V1_post_points = V1_post_points+1;
        end
        if R_duration >= 50*d_adj || R_peak > 1*a_adj
            V1_post_points = V1_post_points+2;

        elseif R_duration >= 40*d_adj || R_peak >= 0.7*a_adj
            V1_post_points = V1_post_points+1;
        end
    
    elseif No_confounders == 1
        if (abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude >= 1)
            V1_post_points = V1_post_points+1;
        end
        if R_duration >= 50*d_adj || R_peak > 1*a_adj
            V1_post_points = V1_post_points+2;

        elseif R_duration >= 40*d_adj || R_peak >= 0.7*a_adj
            V1_post_points = V1_post_points+1;
        end
        
    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.V1, QRST_pts, 5);
        if (abs(Sprime_peak) > 1e-6) && (S_peak/Sprime_peak >= 2.0)
            V1_post_points = V1_post_points + 3;

        elseif (abs(Sprime_peak) > 1e-6) && (S_peak/Sprime_peak >= 1.5)
            V1_post_points = V1_post_points + 2;

        elseif (abs(Sprime_peak) > 1e-6) && (S_peak/Sprime_peak >= 1.25)
            V1_post_points = V1_post_points + 1;
        end
    else
        V1_post_points = 0;
    end

    fprintf('checkpoint here Lead V1 Post\n');
    %Lead V2 Ant
    V2_ant_points = 0;
    [Q_exist] = Q_wave_detection(ecg.V2, QRST_pts);
    
    if Q_exist ~= 0
        
        [Q_duration, Q_peak] = extractQDuration(ecg.V2, QRST_pts, 5);
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.V2, QRST_pts, Q_duration, RBBB_status);
    
    if R_exist && Q_peak
        R_Q = R_peak/Q_peak;

    else
        R_Q = 0;
    end

    [R_duration, initR_duration, initR_peak_amp, initR_peak_loc] = detectRWave(ecg.V2, QRST_pts, Q_duration, 5);
    [hasNtchInit40] = detectNtchInit40(ecg.V2, QRST_pts, 5);

    if RBBB_status == 1 && LAFB_status == 0
        if Q_duration >= 50*d_adj
            V2_ant_points = V2_ant_points +2;
        end
        if Q_exist == 1 || R_duration <= 10*d_adj || R_peak <= 0.1*a_adj
            V2_ant_points = V2_ant_points + 1;
        end

    elseif LAFB_status == 1 && RBBB_status == 0
        
        if (Q_exist == 1 && R_exist == 1) || R_duration <= 10*d_adj || R_peak <= 0.1*a_adj
            V2_ant_points = V2_ant_points + 1;
        end

    elseif LAFB_status == 1 && RBBB_status == 1
        if Q_duration >= 50*d_adj
            V2_ant_points = V2_ant_points + 2;
        end

        if Q_exist == 1 || R_duration <= 10*d_adj || R_peak >= 0.1*a_adj
            V2_ant_points = V2_ant_points + 1;
        end

    elseif LVH_status == 1
        if Q_exist == 1 && R_exist == 1 || hasNtchInit40
            V2_ant_points = V2_ant_points + 1;
        end

    elseif No_confounders == 1
        if Q_exist || R_duration <= 10*d_adj || R_peak <= 0.1*a_adj
            V2_ant_points = V2_ant_points + 1;
        end

    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.V2, QRST_pts, 5);
        if hasNchInit40 == 1
            V2_ant_points = V2_ant_points + 1;
        end
        if R_duration>= 30*d_adj || R_peak >= 0.4*a_adj
            V2_ant_points = V2_ant_points + 2;

        elseif R_duration>= 20*d_adj || R_peak >= 0.3*a_adj
            V2_ant_points = V2_ant_points + 1;
    
        end

    else
        V2_ant_points = 0;
    end 
    
    
    fprintf('checkpoint here Lead V2 Ant\n');
    %Lead V2 post
    V2_post_points = 0;
    [S_amplitude] = detectSWave(ecg.V2, QRST_pts, initR_peak_loc);

    if RAO_status == 1
        a=0;
    elseif RBBB_status == 1 && LAFB_status == 0
        if initR_duration>= 70*d_adj || initR_peak_amp >= 2.5*a_adj
            V2_post_points = V2_post_points + 2;

        elseif initR_duration >= 50*d_adj || initR_peak_amp >= 2.0*a_adj
            V2_post_points = V2_post_points +1;
        end

    elseif LAFB_status == 1 && RBBB_status == 0
        if (abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude >= 1)
            V2_post_points = V2_post_points+1;
        end

        if R_duration >= 60*d_adj || R_peak >= 2*a_adj
            V2_post_points = V2_post_points + 2;

        elseif R_duration >= 50*d_adj || R_peak >= 1.5*a_adj
            V2_post_points = V2_post_points+1;
        end

        if Q_peak <= 0.3*a_adj && S_amplitude <= 0.3*a_adj
            V2_post_points = V2_post_points +1;
        end

    elseif LAFB_status == 1 && RBBB_status == 1
        if initR_duration >= 70*d_adj || initR_peak_amp >= 2.5*a_adj
            V2_post_points = V2_post_points+2;

        elseif initR_duration >= 50*d_adj || initR_peak_amp >= 2.0*a_adj
            V2_post_points = V2_post_points+1;
        end

    elseif LVH_status == 1
        if (abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude >= 1)
            V2_post_points = V2_post_points+1;
        end

        if R_duration >= 60*d_adj || R_peak >= 2*a_adj
            V2_post_points = V2_post_points+2;

        elseif R_duration >= 50*d_adj || R_peak >= 1.5*a_adj
            V2_post_points = V2_post_points+1;
        end
        if Q_peak <= 0.3*a_adj && S_amplitude <= 0.3*a_adj
            V2_post_points = V2_post_points+1;
        end

    elseif No_confounders == 1
        if (abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude >= 1)
            V2_post_points = V2_post_points+1;
        end

        if R_duration >= 60*d_adj || R_peak >= 2*a_adj
            V2_post_points = V2_post_points+2;

        elseif R_duration >= 50*d_adj || R_peak >= 1.5*a_adj
            V2_post_points = V2_post_points+1;
        end

        if Q_peak <= 0.3*a_adj && S_amplitude <= 0.3*a_adj
            V2_post_points = V2_post_points+1;
        end

    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.V2, QRST_pts, 5);
        if (abs(Sprime_peak) > 1e-6) && (S_peak/Sprime_peak >= 2.5) 
            V2_post_points = V2_post_points + 3;

        elseif (abs(Sprime_peak) > 1e-6) && (S_peak/Sprime_peak >= 2.0)
            V2_post_points = V2_post_points + 2;

        elseif (abs(Sprime_peak) > 1e-6) && (S_peak/Sprime_peak >= 1.5)
            V2_post_points = V2_post_points + 1;
        end
    else
        V2_post_points = 0;
    end

    fprintf('checkpoint here Lead V2 post\n');
    %Lead V3
    V3_points = 0;
    [Q_exist] = Q_wave_detection(ecg.V3, QRST_pts);
    
    if Q_exist ~= 0
        
        [Q_duration, Q_peak] = extractQDuration(ecg.V3, QRST_pts, 5);
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.V3, QRST_pts, Q_duration, RBBB_status);
    
    if R_exist && Q_peak
        R_Q = R_peak/Q_peak;

    else
        R_Q = 0;
    end

    [R_duration, initR_duration, initR_peak_amp, initR_peak_loc] = detectRWave(ecg.V3, QRST_pts, Q_duration, 5);
    [hasNtchInit40] = detectNtchInit40(ecg.V3, QRST_pts, 5);

    if RBBB_status == 1 && LAFB_status == 0
        if Q_duration >= 30*d_adj || R_duration <= 10*d_adj
            V3_points = V3_points + 2;

        elseif Q_duration >= 20*d_adj || R_duration <= 20*d_adj
            V3_points = V3_points + 1;
        end

    elseif LAFB_status == 1 && RBBB_status == 0
        if Q_duration >= 30*d_adj || R_duration <= 10*d_adj
            V3_points = V3_points + 2;

        elseif Q_duration >= 20*d_adj || R_duration <= 20*d_adj
            V3_points = V3_points + 1;
        end
    elseif LAFB_status == 1 && RBBB_status == 1
        if Q_duration >= 30*d_adj || R_duration <= 10*d_adj
            V3_points = V3_points + 2;

        elseif Q_duration >= 20*d_adj || R_duration <= 20*d_adj
            V3_points = V3_points + 1;
        end

    elseif LVH_status == 1
        if Q_exist == 1 && R_exist == 1 && Q_duration >= 30*d_adj
            V3_points = V3_points + 2;
        elseif hasNtchInit40 == 1 || (Q_exist == 1 && R_exist == 1) || Q_exist == 1
            V3_points = V3_points + 1;
        end
    elseif No_confounders == 1
        if Q_duration >= 30*d_adj || R_duration <= 10*d_adj
            V3_points = V3_points + 2;

        elseif Q_duration >= 20*d_adj || R_duration <= 20*d_adj
            V3_points = V3_points + 1;
        end
    end
    
    fprintf('checkpoint here Lead V3\n');
    %Lead V4
    V4_points = 0;
    [Q_exist] = Q_wave_detection(ecg.V4, QRST_pts);
    
    if Q_exist ~= 0
        
        [Q_duration, Q_peak] = extractQDuration(ecg.V4, QRST_pts, 5);
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.V4, QRST_pts, Q_duration, RBBB_status);
    
    if R_exist && Q_peak
        R_Q = R_peak/Q_peak;

    else
        R_Q = 0;
    end
    
    [R_duration, initR_duration, initR_peak_amp, initR_peak_loc] = detectRWave(ecg.V4, QRST_pts, Q_duration, 5);
    [S_amplitude] = detectSWave(ecg.V4, QRST_pts, initR_peak_loc);
    [hasNtchInit40] = detectNtchInit40(ecg.V4, QRST_pts, 5);
    
    if RBBB_status == 1 && LAFB_status == 0
        
        if Q_duration >= 20*d_adj
            V4_points = V4_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 0.5) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 0.5)
            V4_points = V4_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1) || R_peak <= 0.5*a_adj || hasNtchInit40 == 1
            V4_points = V4_points + 1;
        end
        
    elseif LAFB_status == 1 && RBBB_status == 0
        
        if Q_duration >= 20*d_adj
            V4_points = V4_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 0.5) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 0.5)
            V4_points = V4_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1) || R_peak <= 0.5*a_adj || hasNtchInit40 == 1
            V4_points = V4_points + 1;
        end
    elseif LAFB_status == 1 && RBBB_status == 1
        
        if Q_duration >= 20*d_adj
            V4_points = V4_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 0.5) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 0.5)
            V4_points = V4_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1) || R_peak <= 0.5*a_adj || hasNtchInit40 == 1
            V4_points = V4_points + 1;
        end
    elseif LVH_status == 1
        
        if Q_duration >= 20*d_adj
            V4_points = V4_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 0.5) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 0.5)
            V4_points = V4_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1) || R_peak <= 0.5*a_adj || hasNtchInit40 == 1
            V4_points = V4_points + 1;
        end
    elseif No_confounders == 1
        
        if Q_duration >= 20*d_adj
            V4_points = V4_points + 1;
        end
        
        if (Q_peak ~= 0 && R_peak/Q_peak <= 0.5) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 0.5)
            
            V4_points = V4_points + 2;
        
        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1) || R_peak <= 0.5*a_adj || hasNtchInit40 == 1
            
            V4_points = V4_points + 1;
        end
    end

    fprintf('checkpoint here Lead V4\n');
    %Lead V5
    V5_points = 0;
    [Q_exist] = Q_wave_detection(ecg.V5, QRST_pts);
    
    if Q_exist ~= 0
        
        [Q_duration, Q_peak] = extractQDuration(ecg.V5, QRST_pts, 5);
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.V5, QRST_pts, Q_duration, RBBB_status);
    
    if R_exist && Q_peak
        R_Q = R_peak/Q_peak;

    else
        R_Q = 0;
    end
    
    [R_duration, initR_duration, initR_peak_amp, initR_peak_loc] = detectRWave(ecg.V5, QRST_pts, Q_duration, 5);
    
    [S_amplitude] = detectSWave(ecg.V5, QRST_pts, initR_peak_loc);
    
    [hasNtchInit40] = detectNtchInit40(ecg.V5, QRST_pts, 5);
    
    if RBBB_status == 1 && LAFB_status == 0
        if Q_duration >= 30*d_adj
            V5_points = V5_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1)
            V5_points = V5_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 2) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 2) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V5_points = V5_points + 1;
        end
    elseif LAFB_status == 1 && RBBB_status == 0
        if Q_duration >= 30*d_adj
            V5_points = V5_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1)
            V5_points = V5_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 2) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1.5) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V5_points = V5_points + 1;
        end
    elseif LAFB_status == 1 && RBBB_status == 1
        if Q_duration >= 30*d_adj
            V5_points = V5_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1)
            V5_points = V5_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 2) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1.5) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V5_points = V5_points + 1;
        end
    elseif LVH_status == 1
        if Q_duration >= 30*d_adj
            V5_points = V5_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1)
            V5_points = V5_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 2) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 2) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V5_points = V5_points + 1;
        end
    elseif No_confounders == 1
        if Q_duration >= 30*d_adj
            V5_points = V5_points + 1;
        end
        if (Q_peak ~= 0 && R_peak/Q_peak <= 1) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 1)
            V5_points = V5_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 2) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 2) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V5_points = V5_points + 1;
        end

    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.V5, QRST_pts, 5);
        if ~Q_peak 
            V5_points = V5_points + 1;
        end
        if (abs(Rprime_peak) > 1e-6) && (R_peak/Rprime_peak >= 2)
            V5_points = V5_points + 2;

        elseif ((abs(Rprime_peak) > 1e-6) && (R_peak/Rprime_peak >= 1)) || ((abs(S_peak) > 1e-6) && (R_peak/S_peak <= 2))
            V5_points = V5_points + 1;
        end

        if R_peak <= 0.5*a_adj
            V5_points = V5_points + 1;
        end
    else
        V5_points = 0;
    end

    fprintf('checkpoint here Lead V5\n');
    %Lead V6
    V6_points = 0;
    [Q_exist] = Q_wave_detection(ecg.V6, QRST_pts);
    
    if Q_exist ~= 0

        [Q_duration, Q_peak] = extractQDuration(ecg.V6, QRST_pts, 5);
        
    else
        Q_duration = 0;
        Q_peak = 0;
    end
    
    [R_exist, R_peak] = R_wave_detection(ecg.V6, QRST_pts, Q_duration, RBBB_status);
    
    if R_exist && Q_peak
        R_Q = R_peak/Q_peak;
        
    else
        R_Q = 0;
    end
    
    [R_duration, initR_duration, initR_peak_amp, initR_peak_loc] = detectRWave(ecg.V6, QRST_pts, Q_duration, 5);
    [S_amplitude] = detectSWave(ecg.V6, QRST_pts, initR_peak_loc);
    [hasNtchInit40] = detectNtchInit40(ecg.V6, QRST_pts, 5);
    
    if RBBB_status == 1 && LAFB_status == 0
        if Q_duration >= 30*d_adj
            V6_points = V6_points + 1;
        end
        if ((abs(Q_peak) > 1e-6) && (R_peak/Q_peak <= 1)) || ((abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude <= 1))
            V6_points = V6_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 3) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 3) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V6_points = V6_points + 1;
        end 
        
    elseif LAFB_status == 1 && RBBB_status == 0
        if Q_duration >= 30*d_adj
            V6_points = V6_points + 1;
        end
        if ((abs(Q_peak) > 1e-6) && (R_peak/Q_peak <= 1)) || ((abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude <= 1))
            V6_points = V6_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 3) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 2) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V6_points = V6_points + 1;
        end 
    elseif LAFB_status == 1 && RBBB_status == 1
        if Q_duration >= 30*d_adj
            V6_points = V6_points + 1;
        end
        if ((abs(Q_peak) > 1e-6) && (R_peak/Q_peak <= 1)) || ((abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude <= 1))
            V6_points = V6_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 3) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 2) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V6_points = V6_points + 1;
        end 
    elseif LVH_status == 1
        if Q_duration >= 30*d_adj
            V6_points = V6_points + 1;
        end
        if ((abs(Q_peak) > 1e-6) && (R_peak/Q_peak <= 1)) || ((abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude <= 1))
            V6_points = V6_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 3) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 3) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V6_points = V6_points + 1;
        end 
    elseif No_confounders == 1
        if Q_duration >= 30*d_adj
            V6_points = V6_points + 1;
        end
        if ((abs(Q_peak) > 1e-6) && (R_peak/Q_peak <= 1)) || ((abs(S_amplitude) > 1e-6) && (R_peak/S_amplitude <= 1))
            V6_points = V6_points + 2;

        elseif (Q_peak ~= 0 && R_peak/Q_peak <= 3) || (S_amplitude ~= 0 && R_peak/S_amplitude <= 3) || R_peak <= 0.6*a_adj || hasNtchInit40 == 1
            V6_points = V6_points + 1;
        end 
    elseif LBBB_status == 1
        [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(ecg.V6, QRST_pts, 5);
        if Q_duration >= 20*d_adj
            V6_points = V6_points + 1;
        end
        if (abs(Rprime_peak) > 1e-6) && (R_peak/Rprime_peak >= 2)
            V6_points = V6_points + 2;

        elseif ((abs(Rprime_peak) > 1e-6) && (R_peak/Rprime_peak >= 1)) || ((abs(S_peak) > 1e-6) && (R_peak/S_peak <= 2))
            V6_points = V6_points + 1;
        end

        if R_peak <= 0.6*a_adj
            V6_points = V6_points + 1;
        end
    else
        V6_points = 0;
    end
    
    fprintf('checkpoint here Lead V6\n');
    score = L1_points+L2_points+avL_points+avF_points+V1_ant_points+V1_post_points+V2_ant_points+V2_post_points+...
        V3_points+V4_points+V5_points+V6_points;

    end


function [Q_peak, Q_duration, R_peak, R_duration, Rprime_peak, S_peak, Sprime_peak, hasNchInit40] = getLBBBMeasurements(signal, QRST_pts, window_size)
    % Initialize outputs
    Q_peak = 0; Q_duration = 0;
    R_peak = 0; R_duration = 0; Rprime_peak = 0;
    S_peak = 0; Sprime_peak = 0;
    hasNchInit40 = 0;
    
    QRS_start = QRST_pts(1);
    QRS_end = QRST_pts(3);
    baseline_window = 20;
    baseline = median(signal(QRST_pts(5):min(QRST_pts(5) + baseline_window, length(signal))));
    ms40_point = QRS_start + 40;  % Assuming 1000Hz
    
    % Q wave detection
    if signal(QRS_start) > signal(QRS_start+5)  % Check negative deflection
        [Q_peak, ~] = min(signal(QRS_start:QRST_pts(2)));
        % Find Q duration (until baseline crossing)
        for i = QRS_start+5:QRST_pts(2)
            if signal(i) >= baseline
                Q_duration = i - QRS_start;
                break;
            end
        end
    end
    
    % R wave and R' wave detection
    dv1 = diff(signal(QRS_start+Q_duration:QRS_end));
    dv1_smooth = movmean(dv1, window_size);
    r_peaks = [];

    % Find all positive peaks
    for i = window_size+1:length(dv1_smooth)-window_size
        % Check for positive to negative crossing with window confirmation
        if mean(dv1_smooth(i-window_size+1:i)) > 0 && ...  % Previous window positive
           mean(dv1_smooth(i+1:i+window_size)) < 0         % Next window negative
            peak_loc = QRS_start + i;
            peak_amp = signal(peak_loc);
            if peak_amp > baseline  % Only count peaks above baseline
                r_peaks = [r_peaks; peak_loc peak_amp];
            end
        end
    end

    % First peak is R, next significant peak is R'
    if ~isempty(r_peaks)
        R_peak = r_peaks(1,2);  % Amplitude of first peak
        r_loc = r_peaks(1,1);   % Location of first peak
        
        % Calculate R duration (from QRS start to baseline crossing)
        for i = r_loc:QRS_end
            if signal(i) <= baseline
                R_duration = i - QRS_start;
                break;
            end
        end
        
        % Look for next significant peak after 40ms for R'
        for i = 2:size(r_peaks,1)
            if r_peaks(i,1) > ms40_point
                Rprime_peak = r_peaks(i,2);
                break;
            end
        end
    end
    
    % S wave and S' wave detection
    s_peaks = [];
    
    % Find all negative peaks after R wave
    for i = (r_loc - QRS_start + window_size):length(dv1_smooth)-window_size
        % Check for negative to positive crossing with window confirmation
        if mean(dv1_smooth(i-window_size+1:i)) < 0 && ...  % Previous window negative
           mean(dv1_smooth(i+1:i+window_size)) > 0         % Next window positive
            peak_loc = QRS_start + i;
            peak_amp = signal(peak_loc);
            if peak_amp < baseline  % Only count troughs below baseline
                s_peaks = [s_peaks; peak_loc peak_amp];
            end
        end
    end
    
    % First trough is S, next significant trough is S'
    if ~isempty(s_peaks)
        S_peak = s_peaks(1,2);  % Amplitude of first trough
        % Look for next significant trough after 40ms for S'
        for i = 2:size(s_peaks,1)
            if s_peaks(i,1) > ms40_point
                Sprime_peak = s_peaks(i,2);
                break;
            end
        end
    end
    
    % NchInit40 detection (in first 40ms)
    dv1_40ms = diff(signal(QRS_start:ms40_point));
    dv1_smooth_40ms = movmean(dv1_40ms, window_size);
    
    for i = window_size+1:length(dv1_smooth_40ms)-window_size
        v1 = dv1_smooth_40ms(i-window_size:i);
        v2 = dv1_smooth_40ms(i:i+window_size);
        
        vec1 = mean(v1);
        vec2 = mean(v2);
        angle_deg = acosd(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));
        
        if angle_deg >= 90
            hasNchInit40 = 1;
            break;
        end
    end
end


function [S_amplitude] = detectSWave(signal, QRST_pts, initR_end)
    % Initialize output
    S_amplitude = 0;
    
    % Define window from initial R end to QRS end
    QRS_end = QRST_pts(3);
    baseline_window = 20;
    baseline = median(signal(QRST_pts(5):min(QRST_pts(5) + baseline_window, length(signal))));
    
    % Find where signal crosses baseline after initial R
    S_start = initR_end;
    if initR_end == 0
        S_start = S_start+1;

    end

    for i = S_start:QRS_end
        if signal(i) <= baseline
            S_start = i;
            break;
        end
    end

    % Return if no valid window
    if S_start >= QRS_end
        return;
    end
    
    % Find S wave amplitude (minimum after crossing baseline)
    [S_amplitude, ~] = min(signal(S_start:QRS_end));
end

function [hasNtchInit40] = detectNtchInit40(signal, QRST_pts, window_size)
    % Use existing angle detection code, but only look at first 40ms
    QRS_start = QRST_pts(1);
    ms40_idx = QRS_start + 40;  % Assuming 1000Hz sampling rate, adjust if different
    
    % Use our existing derivative and angle calculation
    dv1 = diff(signal(QRS_start:ms40_idx));
    dv1_smooth = movmean(dv1, window_size);
    
    hasNtchInit40 = 0;
    
    % Look for direction changes in first 40ms only
    for i = window_size+1:length(dv1_smooth)-window_size
        v1 = dv1_smooth(i-window_size:i);
        v2 = dv1_smooth(i:i+window_size);
        
        vec1 = mean(v1);
        vec2 = mean(v2);
        angle_deg = acosd(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));
        
        if angle_deg >= 90
            hasNtchInit40 = 1;
            break;
        end
    end
end


function [R_duration, initR_duration, initR_peak_amp, initR_peak_loc] = detectRWave(signal, QRST_pts, Q_duration, window_size)
    % Initialize outputs
    R_duration = 0;
    initR_duration = 0;
    initR_peak_amp = 0;
    initR_peak_loc = 0;
    
    % Get QRS window
    QRS_start = QRST_pts(1);
    QRS_end = QRST_pts(3);

    baseline_window = 20;
    baseline = median(signal(QRST_pts(5):min(QRST_pts(5) + baseline_window, length(signal))));
    
    % If Q wave present, start from where Q returns to baseline
    if Q_duration > 0
        R_start = QRS_start + Q_duration;
        % Find where Q returns to baseline
        for i = R_start:QRS_end
            if signal(i) >= baseline
                R_start = i;
                break;
            end
        end
    else
        R_start = QRS_start;  % No Q wave
    end
    
    % Calculate derivatives for R wave analysis
    dv1 = diff(signal(R_start:QRS_end));
    dv1_smooth = movmean(dv1, window_size);
    
    % Find first positive deflection
    found_R = 0;
    for i = window_size:length(dv1_smooth)-window_size
        % Check if average slope over window is positive
        window_avg = mean(dv1_smooth(i-window_size+1:i+window_size));
        if window_avg > 0
            found_R = 1;
            R_start = R_start + i;
            break;
        end
    end
    
    if ~found_R
        return;  % No R wave found
    end
    
    % Look for notch (for initR)
    hasNotch = 0;
    notch_peak = 0;
    
    for i = window_size+1:length(dv1_smooth)-window_size
        v1 = dv1_smooth(i-window_size:i);
        v2 = dv1_smooth(i:i+window_size);
        
        vec1 = mean(v1);
        vec2 = mean(v2);
        angle_deg = acosd(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));
        
        if angle_deg >= 90
            hasNotch = 1;
            notch_peak = R_start + i;
            break;
        end
    end
    
    % Find R peak and calculate durations
    if hasNotch
        % For initR calculations
        [initR_peak_amp, peak_idx] = max(signal(R_start:notch_peak));
        initR_peak_loc = R_start + peak_idx - 1;
        
        % Check 50ms requirement for initR
        if (initR_peak_loc - QRS_start) <= 50
            initR_duration = notch_peak - R_start;
        end
        
        % For regular R duration, find where signal crosses baseline after notch
        for i = notch_peak:QRS_end
            if signal(i) <= baseline
                R_duration = i - R_start;
                break;
            end
        end
    else
        % Find where signal returns to baseline
        [peak_amp, peak_idx] = max(signal(R_start:QRS_end));
        peak_loc = R_start + peak_idx - 1;
        
        % Check 50ms for initR
        if (peak_loc - QRS_start) <= 50
            initR_peak_amp = peak_amp;
            initR_peak_loc = peak_loc;
        end
        
        % Find end of R wave at baseline crossing
        for i = peak_loc:QRS_end
            if signal(i) <= baseline
                R_duration = i - R_start;
                initR_duration = R_duration;  % Same for non-notched R
                break;
            end
        end
    end
end



function [R_status, R_peak] = R_wave_detection(signal, QRST_pts, Q_duration, RBBB_status)
    R_status = 0;
    R_peak = 0;

    R_start = QRST_pts(1) + Q_duration;
    R_end = QRST_pts(3);

    % Calculate derivative for search window
    search_signal = signal(R_start:R_end);
    dv1 = diff(search_signal);
    window_size = 5;
    dv1_smooth = movmean(dv1, window_size);
    
    % Look for first positive deflection in window
    for i = window_size:length(dv1_smooth)-window_size
        if dv1_smooth(i) > 0  % Found positive slope
            % Track until slope changes to negative to find peak
            for j = i:length(dv1_smooth)
                if dv1_smooth(j) < 0  % Found peak
                    R_status = 1;
                    return;
                end
            end
            break;  % Only look for first positive deflection
        end
    end

    if RBBB_status == 0
        R_peak = max(search_signal);

        baseline_window = 20;
        baseline = median(signal(QRST_pts(5):min(QRST_pts(5) + baseline_window, length(signal))));
        R_peak = R_peak - baseline;
    end
    
end

function [Q_status] = Q_wave_detection(signal, QRST_pts)
    if signal(QRST_pts(1)) >= signal(QRST_pts(1) + 10)
        %go 10 samples to the right and see
        Q_status = 1;
    else
        Q_status = 0;
    
    end
end


function [Q_duration, Q_peak] = extractQDuration(signal, QRST_pts, window_size)
    
    Q_duration = 0;
    Q_peak = 0;
    % Calculate first derivative
    qrs_signal = signal(QRST_pts(1):QRST_pts(3));
    dv1 = diff(qrs_signal);
    dv1_smooth = movmean(dv1, window_size);
    
    % Use QRST onset as Q start
    Q_start = QRST_pts(1);
    
    % Initialize variables relative to Q_start
    Q_end = Q_start;
    hasNotch = 0;
    notch_peak = 0;
    
    % Look for direction changes indicating notch
    % Adjust loop to use global indices

    baseline_window = 20;
    baseline = median(signal(QRST_pts(5):min(QRST_pts(5) + baseline_window, length(signal))));
    
    for i = 1+window_size:length(dv1_smooth)-window_size
        % Get vectors before and after point, using global indices
        v1 = dv1_smooth(i-window_size:i);
        v2 = dv1_smooth(i:i+window_size);
        
        % Rest of notch detection remains same
        vec1 = mean(v1);
        vec2 = mean(v2);
        angle_deg = acosd(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));
        
        if angle_deg >= 90 && abs(qrs_signal(i+1) - qrs_signal(i)) >= 0.05
            hasNotch = 1;
            notch_peak = Q_start+i;  % This will now be a global index
            [Q_peak, ~] = min(signal(Q_start:notch_peak));
            break;
        end
    end
    
    if hasNotch
        Q_end = notch_peak;
        
    else
        
        for i = QRST_pts(1):QRST_pts(3)
            if signal(i) >= baseline && i-Q_start >=5
                Q_end = i;
                break;
            end
        end
        
        [Q_peak, ~] = min(signal(Q_start:Q_end));
    end
    
    Q_duration = Q_end - Q_start;
end


function [RAO_status] = check_RAO(ecg, QRST_pts)
    V1_signal = ecg.V1;
    avF_signal = ecg.avF;
    
    QRS_onset = QRST_pts(1);
    %We dont technically need 
    pre_qrs_max_V1 = max(V1_signal(1:QRS_onset));
    pre_qrs_max_avF = max(avF_signal(1:QRS_onset));
    
    if pre_qrs_max_V1 >= 0.1 || pre_qrs_max_avF >= 0.175
        RAO_status = 1;
    else
        RAO_status = 0;
    end
end

function [LVH_status] = check_LVH(ecg, isMale, fidpts)
    %Lead V3
    V3_signal = ecg.V3;
    [~, V3_s_wave, ~, ~, ~] = rs_values(V3_signal,fidpts);

    %Lead avL
    avL_signal = ecg.avL;
    [avL_r_wave, ~, ~, ~, ~] = rs_values(avL_signal,fidpts);

    %Lead V1
    V1_signal = ecg.V1;
    [~, V1_s_wave, ~, ~, ~] = rs_values(V1_signal,fidpts);

    %Lead V5
    V5_signal = ecg.V5;
    [V5_r_wave, ~, ~, ~, ~] = rs_values(V5_signal,fidpts);

    %Lead V6
    V6_signal = ecg.V6;
    [V6_r_wave, ~, ~, ~, ~] = rs_values(V6_signal,fidpts);
    
    cornell_lvh_mv = abs(V3_s_wave) + avL_r_wave;

    sokolow_lvh_mv_1 = abs(V1_s_wave) + max(V5_r_wave, V6_r_wave);
    sokolow_lvh_mv_2 = max(V5_r_wave, V6_r_wave);
    
    if sokolow_lvh_mv_1 >= 3.50
        LVH_status = 1;
    elseif cornell_lvh_mv >= 2.80 && isMale && sokolow_lvh_mv_2 > 2.60
        LVH_status = 1;
    elseif cornell_lvh_mv >= 2.00 && ~isMale && sokolow_lvh_mv_2 > 2.60
        LVH_status = 1;
    else
        LVH_status = 0;
    end
        
end



function [LAFB_status] = check_LAFB(ecg, QRST_pts, isMale, fidpts)
    QRS_duration = QRST_pts(3) - QRST_pts(1);

    if isMale && QRS_duration >= 100
        QRS_criterion = 1;
    elseif ~isMale && QRS_duration >= 90
        QRS_criterion = 1;
    else
        QRS_criterion = 0;
    end
    
    %Lead 1
    L1_signal = ecg.I;
    [L1_r_wave, L1_s_wave, ~, ~, ~] = rs_values(L1_signal,fidpts);
    
    %Lead avF
    avF_signal = ecg.avF;
    [avF_r_wave, avF_s_wave, ~, ~, ~] = rs_values(avF_signal,fidpts);

    qrs_frontal_axis = ecg_axis(L1_r_wave, L1_s_wave, avF_r_wave, avF_s_wave);

    if qrs_frontal_axis <= -45 && qrs_frontal_axis > -180
        left_axis = 1;
    else
        left_axis = 0;
    end

    if left_axis == 1 && QRS_criterion == 1
        LAFB_status = 1;
    else
        LAFB_status = 0;
    end

end


function [terminal_deflection] = analyzeV1WavePattern(ECG, QRST_pts)
    % Get V1 signal and QRS complex
    v1_signal = ECG.V1;
    waveType = "None";
    qrs_complex = v1_signal(QRST_pts(1):QRST_pts(3));
    
    % Get baseline from segment after T wave end
    % Use 20ms (about 20 samples at 1000Hz) after T end
    baseline_window = 20;
    baseline = median(v1_signal(QRST_pts(5):min(QRST_pts(5) + baseline_window, length(v1_signal))));
    
    % Define terminal portion (last 40% of QRS)
    terminal_percent = 0.4;
    terminal_start = floor(length(qrs_complex) * (1-terminal_percent)) + 1;
    terminal_portion = qrs_complex(terminal_start:end);
    
    % Calculate and smooth derivative of terminal portion
    dv1 = diff(terminal_portion);
    window_size = 5;  % 5ms at 1000Hz
    dv1_smooth = movmean(dv1, window_size);
    
    % Find direction changes in terminal portion
    direction_changes = [];
    i = window_size;
    while i <= length(dv1_smooth)-window_size
        if mean(dv1_smooth(i-window_size+1:i)) > 0 && ...
           mean(dv1_smooth(i+1:i+window_size)) < 0
            direction_changes = [direction_changes; i];
            i = i + window_size;
        elseif mean(dv1_smooth(i-window_size+1:i)) < 0 && ...
               mean(dv1_smooth(i+1:i+window_size)) > 0
            direction_changes = [direction_changes; i];
            i = i + window_size;
        else
            i = i + 1;
        end
    end
    
    % Find peaks and troughs in terminal portion
    peaks = [];
    troughs = [];
    for i = 1:length(direction_changes)
        if i == 1
            prev_window = 1:direction_changes(i);
        else
            prev_window = direction_changes(i-1):direction_changes(i);
        end
        if mean(dv1_smooth(prev_window)) > 0
            peaks = [peaks; direction_changes(i)];
        else
            troughs = [troughs; direction_changes(i)];
        end
    end
    
    % Get amplitudes relative to baseline
    peak_amplitudes = terminal_portion(peaks) - baseline;
    trough_amplitudes = terminal_portion(troughs) - baseline;
    
    % Analyze terminal wave pattern
    if length(peaks) == 1 && isempty(troughs)
        
        if peak_amplitudes(1) > 0
            waveType = 'R';
        end
    elseif isempty(peaks) && length(troughs) == 1
        
        if trough_amplitudes(1) < 0
            waveType = 'Q';
        end 

    elseif ~isempty(peaks) && length(peaks) > 1
        
        waveType = 'R-prime';
    else
        
        if ~isempty(peaks) && ~isempty(troughs)
            first_peak_idx = peaks(1);
            first_trough_idx = troughs(1);
            
            if first_peak_idx < first_trough_idx && ...
               abs(peak_amplitudes(1)) < abs(trough_amplitudes(1))
                waveType = 'rS';
            else
                waveType = 'R';
            end
        end
        
    end
    
    if waveType == "Q" || waveType == "rS"
        
        terminal_deflection = 0;
    elseif waveType == "R" || waveType == "R-prime"
        
        terminal_deflection = 1;
    else
        
        terminal_deflection = 1;
    end

end

function [LBBB_status] = check_LBBB(ecg, QRST, isMale)
    QRS_durations = QRST(3) - QRST(1);
    if isMale
        QRS_criterion = QRS_durations >=140;
    else
        QRS_criterion = QRS_durations >=130;
    end
    
    mid_qrs_start = QRST(1) + 40;
    leads_to_check = {"V1", "V2", "V5", "V6", "I", "avL"};
    leads_with_pattern = {};
    window_size = 5;
    
    for i=1:length(leads_to_check)
        lead = leads_to_check{i};
        signal = ecg.(lead);
        
        % Check if indices are valid
        if mid_qrs_start > length(signal) || QRST(3) > length(signal) || mid_qrs_start >= QRST(3)
            fprintf('Warning: Invalid indices for lead %s\n', lead);
            continue;
        end
        
        mid_qrs = signal(mid_qrs_start:QRST(3));
        
        % Check for slurring and notching
        dv1 = diff(mid_qrs);
        dv1_smooth = movmean(dv1, window_size);
        mean_velocity = mean(abs(dv1_smooth));
        low_velocity_threshold = mean_velocity * 0.3;
        
        low_velocity_segments = abs(dv1_smooth) < mean_velocity & ...
                               abs(dv1_smooth) > low_velocity_threshold;
        
        min_duration = 5;
        
        % Fix for has_slurring calculation
        if length(low_velocity_segments) >= min_duration
            conv_result = conv(double(low_velocity_segments), ones(1, min_duration), 'valid');
            has_slurring = any(conv_result == min_duration);
        else
            has_slurring = false;
        end
        
        % Make sure input to detectQRSNotch is valid
        if length(mid_qrs) > 2*window_size + 1
            [has_notching] = detectQRSNotch(mid_qrs, window_size);
        else
            has_notching = 0;
        end
        
        if has_notching || has_slurring
            leads_with_pattern{end+1} = lead;
        end
    end
    
    if length(leads_with_pattern) >= 2 && QRS_criterion == 1
        LBBB_status = 1;
    else
        LBBB_status = 0;
    end
end

function [hasNotch] = detectQRSNotch(signal, window_size)
    % Calculate first derivative
    dv1 = diff(signal);
    
    % Smooth derivative to reduce noise
    dv1_smooth = movmean(dv1, window_size);
    
    % Initialize outputs
    hasNotch = 0;
    
    % Look for direction changes
    for i = window_size+1:length(dv1_smooth)-window_size
        % Get vectors before and after point
        v1 = dv1_smooth(i-window_size:i); % Vector before point
        v2 = dv1_smooth(i:i+window_size); % Vector after point
        
        % Calculate mean vectors for more stable angle calculation
        vec1 = mean(v1);
        vec2 = mean(v2);
        
        % Avoid division by zero
        if norm(vec1) == 0 || norm(vec2) == 0
            continue;
        end
        
        % Calculate angle between vectors using arctan2
        % angle = arccos(dot(v1,v2)/(norm(v1)*norm(v2)))
        dot_product = dot(vec1, vec2);
        magnitudes_product = norm(vec1) * norm(vec2);
        
        % Ensure dot product is within valid range for acosd
        cos_angle = min(max(dot_product / magnitudes_product, -1), 1);
        angle_deg = acosd(cos_angle);
        
        % Check if angle is ≥90 degrees
        if angle_deg >= 90
            hasNotch = 1;
            break;  % Exit once a notch is found
        end
    end
end