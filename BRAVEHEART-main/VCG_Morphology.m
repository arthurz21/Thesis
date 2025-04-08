%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BRAVEHEART - Open source software for electrocardiographic and vectorcardiographic analysis
% VCG_Morphology.m -- VCG_Morphology Results Class
% Copyright 2016-2024 Hans F. Stabenau and Jonathan W. Waks
% 
% Source code/executables: https://github.com/BIVectors/BRAVEHEART
% Contact: braveheart.ecg@gmail.com
% 
% BRAVEHEART is free software: you can redistribute it and/or modify it under the terms of the GNU 
% General Public License as published by the Free Software Foundation, either version 3 of the License, 
% or (at your option) any later version.
%
% BRAVEHEART is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. 
% If not, see <https://www.gnu.org/licenses/>.
%
% This software is for research purposes only and is not intended to diagnose or treat any disease.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef VCG_Morphology
% Various VCG measurements based on morphology of the VCG loops
    
properties (SetAccess=private)
        
        TCRT                    % Total Cosine R to T
        TCRT_angle              % Angle from TCRT = acos(TCRT)

        tloop_residual          % SVD residual from fitting T loop to a plane (0 = perfect fit) = (t_S3)^2
        tloop_rmse              % RMSE for fit of T loop to best fit plane (0 = perfect fit)
        tloop_roundness         % How round the T loop is.  1 = circle, larger values are elliptical
        tloop_area              % Area of T loop
        tloop_perimeter         % Length of T loop

        qrsloop_residual        % SVD residual from fitting QRS loop to a plane (0 = perfect fit) = (qrs_S3)^2
        qrsloop_rmse            % RMSE for fit of QRS loop to best fit plane (0 = perfect fit)
        qrsloop_roundness       % How round the QRS loop is.  1 = circle, larger values are elliptical
        qrsloop_area            % Area of QRS loop
        qrsloop_perimeter       % Length of QRS loop
      
        
        qrs_S1                  % QRS loop 1st singular value
        qrs_S2                  % QRS loop 2nd singular value
        qrs_S3                  % QRS loop 3rd singular value
        
        t_S1                    % T loop 1st singular value
        t_S2                    % T loop 2nd singular value
        t_S3                    % T loop 3rd singular value
        
        qrs_var_s1_total        % Pct of total variance made up by 1st QRS singular value
        qrs_var_s2_total        % Pct of total variance made up by 2nd QRS singular value
        qrs_var_s3_total        % Pct of total variance made up by 3rd QRS singular value
        
        t_var_s1_total          % Pct of total variance made up by 1st T singular value
        t_var_s2_total          % Pct of total variance made up by 2nd T singular value
        t_var_s3_total          % Pct of total variance made up by 3rd T singular value
        
        qrs_loop_normal         % Vector normal to best fit QRS loop plane
        t_loop_normal           % Vector normal to best fit T loop plane
        qrst_dihedral_ang       % Dihedral angle between best fit QRS loop and T loop planes

        TMD                     % T-Wave Morphology Dispersion (deg)
        TWR_abs                 % Absolute T-wave residuum (mv2)
        TWR_rel                 % Relative T-wave residuum (%)

        qrst_var_s1_total
        qrst_var_s2_total
        qrst_var_s3_total
        qrst_S1
        qrst_S2
        qrst_S3
        qrstloop_roundness
        qrstloop_area
        qrstloop_perimeter
        qrstloop_residual
        qrstloop_rmse
        qrst_loop_normal

        RWR_abs
        RWR_rel

        JWR_abs
        JWR_rel

        STR_abs
        STR_rel

        JTc

        st80_I
        st80_III

        st80_aVL
        st80_aVF

        st80_V1
        st80_V2
        st80_V3
        st80_V4
        st80_V5
        st80_V6

        TMDpre
        TMDpost

        vat_II
        vat_III
        vat_V2
        vat_V4

        jt_anterior
        jt_lateral
        
        latConcaveAmp
        antConcaveAmp

        qdur_avF
        
        true_t_S1
        true_t_S2
        
        TampInfl1
        J_PCAratio

        fpTinfl1Axis
        fpTinfl1Mag

        Tasym
        RMSmin
        RMSvar

        %s1_8lead
        %s2_8lead
        %s3_8lead
        
        MIsz
        PR
end
    
    
methods
        
    function obj = VCG_Morphology(varargin)        % varagin: ECG12, VCG, Beats
        if nargin == 0; return; end
        if nargin ~= 3
            error('VCG_Morphology: expected 3 args in constructor, got %d', nargin);
        end

        assert(isa(varargin{1}, 'ECG12'), 'First argument is not a ECG12 class');
        assert(isa(varargin{2}, 'VCG'), 'Second argument is not a VCG class');
        assert(isa(varargin{3}, 'Beats'), 'Third argument is not a Beats class');

        ecg = varargin{1}; vcg = varargin{2}; fidpts = varargin{3}.beatmatrix();
        
        % TCRT
        [obj.TCRT, obj.TCRT_angle] = tcrt(ecg, vcg, fidpts, 0.7, 0);
        
        % Loop morphology
        [obj.t_loop_normal, ~, ~, ~, obj.tloop_residual, obj.tloop_rmse, obj.t_var_s1_total,  obj.t_var_s2_total,obj.t_var_s3_total, obj.t_S1, obj.t_S2, obj.t_S3, ~, obj.tloop_roundness, ~, ~, ~, obj.tloop_area, obj.tloop_perimeter] = ...
            plane_svd(vcg.X(fidpts(3):fidpts(4)), vcg.Y(fidpts(3):fidpts(4)), vcg.Z(fidpts(3):fidpts(4)), 0);
        
        [obj.qrs_loop_normal, ~, ~, ~, obj.qrsloop_residual, obj.qrsloop_rmse, obj.qrs_var_s1_total,  obj.qrs_var_s2_total,  obj.qrs_var_s3_total, obj.qrs_S1, obj.qrs_S2, obj.qrs_S3, ~, obj.qrsloop_roundness, ~, ~, ~, obj.qrsloop_area, obj.qrsloop_perimeter] = ...
            plane_svd(vcg.X(fidpts(1):fidpts(3)), vcg.Y(fidpts(1):fidpts(3)), vcg.Z(fidpts(1):fidpts(3)), 0);

        [obj.qrst_loop_normal, ~, ~, ~, obj.qrstloop_residual, obj.qrstloop_rmse, obj.qrst_var_s1_total,  obj.qrst_var_s2_total,  obj.qrst_var_s3_total, obj.qrst_S1, obj.qrst_S2, obj.qrst_S3, ~, obj.qrstloop_roundness, ~, ~, ~, obj.qrstloop_area, obj.qrstloop_perimeter] = ...
            plane_svd(vcg.X(fidpts(1):fidpts(4)), vcg.Y(fidpts(1):fidpts(4)), vcg.Z(fidpts(1):fidpts(4)), 0);
        
        % Dihedral angle between QRS and T loop planes
        obj.qrst_dihedral_ang = dihedral(obj.t_loop_normal, obj.qrs_loop_normal);
        
        
        all_points = varargin{3};
        QRST_pts = [all_points.Q all_points.QRS all_points.S all_points.T all_points.Tend];
        
        [R_wave_limits, J_wave_limits, ST_segment_limits] = index_finder_vcg(vcg, QRST_pts);

        [obj.RWR_abs, obj.RWR_rel] = svd_rwave(ecg, R_wave_limits, fidpts);
        
        
        if isnan(J_wave_limits.end) % J_wave is not present
            % obj.JWR_abs = [];
            % obj.JWR_rel = [];
            J_wave_limits.end = J_wave_limits.start+40;
            [obj.JWR_abs, obj.JWR_rel] = svd_jwave(ecg, J_wave_limits);
        else
            [obj.JWR_abs, obj.JWR_rel] = svd_jwave(ecg, J_wave_limits);
        end
        
        [obj.STR_abs, obj.STR_rel] = svd_stwave(ecg, ST_segment_limits);
        
        %[obj.STWR_abs, obj.STWR_rel] = svd_stwave(ST_segment_limits);
        
        % TMD/TWR
        [obj.TMD, obj.TWR_abs, obj.TWR_rel] = svd_twave(ecg,fidpts);
        
        %[obj.ecg] = ecg;
        
        obj.JTc = (fidpts(4) - fidpts(3))*(1000/(ecg.hz));
        
        S_end = fidpts(3);
        st80 = S_end + 80*((ecg.hz)/1000);
        
        obj.st80_I = ecg.I(st80);
        obj.st80_III = ecg.III(st80);

        obj.st80_aVL = ecg.avL(st80);
        obj.st80_aVF = ecg.avF(st80);

        obj.st80_V1 = ecg.V1(st80);
        obj.st80_V2 = ecg.V2(st80);
        obj.st80_V3 = ecg.V3(st80);
        obj.st80_V4 = ecg.V4(st80);
        obj.st80_V5 = ecg.V5(st80);
        obj.st80_V6 = ecg.V6(st80);
        
        t_start = ST_segment_limits.ST_end;
        t_peak = all_points.T;
        t_end = all_points.Tend;
        
        obj.TMDpre = svd_twave_2(ecg, t_start, t_peak);
        
        obj.TMDpost = svd_twave_2(ecg, t_peak, t_end);
        
        [obj.vat_II, obj.vat_III, obj.vat_V2, obj.vat_V4] = find_VAT(ecg, QRST_pts);
        
        [obj.jt_anterior, obj.jt_lateral] = JT_ant_lat_finder(ecg, (ecg.hz)/1000, QRST_pts);
        
        [obj.latConcaveAmp, obj.antConcaveAmp] = computeConcaveAmp(ecg, QRST_pts);
        
        
        [obj.qdur_avF] = qdur_avF(ecg, QRST_pts);

        [~, ~, ~, ~, ~, ~, ~, ~,~, obj.true_t_S1, obj.true_t_S2, ~, ~, ~, ~, ~, ~, ~, ~] = ...
            plane_svd(vcg.X(t_start:fidpts(4)), vcg.Y(t_start:fidpts(4)), vcg.Z(t_start:fidpts(4)), 0);
        
        [obj.TampInfl1, obj.fpTinfl1Axis, obj.fpTinfl1Mag] = calculate_global_first_inflection(vcg, t_start, t_end);
        %fprintf('check point here\n');
        [obj.Tasym] = calculate_asym_score(vcg, t_start, t_end);
        
        [obj.RMSmin, obj.RMSvar] = calculateECGMetrics(ecg, QRST_pts(1),t_end);
        
        
        if ~isnan(J_wave_limits.end)
            [~, ~, ~, ~, ~, ~, ~, ~,~, j_S1, j_S2, ~, ~, ~, ~, ~, ~, ~, ~] = ...
            plane_svd(vcg.X(t_start:fidpts(4)), vcg.Y(t_start:fidpts(4)), vcg.Z(t_start:fidpts(4)), 0);

            obj.J_PCAratio = (j_S2)/(j_S1);
        else
            obj.J_PCAratio = [];
        end
        
        %[s1_8lead, s2_8lead, s3_8lead] = my_own_svd_8lead(ecg)
        %disp('this is ECG data')
        %ecg.I
        
        isMale = 0;
        age = 1;
        [obj.MIsz] = ComputeSzScore(ecg, QRST_pts, fidpts, isMale, age);
        
        %fprintf('check point here\n');
        %[obj.PR] = PR_interval_check_non_PT(ecg, QRST_pts);
        
    end     % End main VCG_Morphology methods 
    
    
%     function values = values(obj)
%     txt_labels = properties(obj);
% 
%     values = zeros(1, length(txt_labels));
%     for i = 1:length(txt_labels)
%         if ~isempty(obj.(txt_labels{i}))
%             values(i) = obj.(txt_labels{i});
%         else
%             values(i) = NaN;
%         end
%     end
%     end
    
function v = cells(obj)
            
            lab = obj.labels();
            v = cell(1, length(lab));
            for i = 1:length(lab)
                
                v{i} = mat2str(obj.(lab{i})');
                T = v{i};
                T(T=='[') = [];
                T(T==']') = [];
                v{i} = T;
                
            end
            
        end    
               
end     % End methods
    
    
methods(Static)
    
    % Display lead header for export purposes
    function labels = labels()
        obj = VCG_Morphology();
        labels = properties(obj)';
    end
    
    % Length of class members
    function l = length(); g = VCG_Morphology(); l = length(properties(g)); end
        
    % Set to all nan if error
    function a = allnan()
        a = VCG_Morphology();
        p = properties(a);
        for i = 1:length(p)
            a.(p{i}) = nan;
        end
    end
   
end     % End methods

end     % End class


