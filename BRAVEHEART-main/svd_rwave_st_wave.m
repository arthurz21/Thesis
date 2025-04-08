
function [R_wave_limits, RWR_abs, RWR_rel] = svd_rwave_st_wave(ecg, fidpts)
%cut out ST segments. we will also be checking for Jwave_here
%we already know S
R_wave_limits = find_R_wave_start_end_with_sign_check(ecg, fidpts);

I_start = round(R_wave_limits.I.R_wave_start);
I_end = round(R_wave_limits.I.R_wave_end);
II_start = round(R_wave_limits.II.R_wave_start);
II_end = round(R_wave_limits.II.R_wave_end);
V1_start = round(R_wave_limits.V1.R_wave_start);
V1_end = round(R_wave_limits.V1.R_wave_end);
V2_start = round(R_wave_limits.V2.R_wave_start);
V2_end = round(R_wave_limits.V2.R_wave_end);
V3_start = round(R_wave_limits.V3.R_wave_start);
V3_end = round(R_wave_limits.V3.R_wave_end);
V4_start = round(R_wave_limits.V4.R_wave_start);
V4_end = round(R_wave_limits.V4.R_wave_end);
V5_start = round(R_wave_limits.V5.R_wave_start);
V5_end = round(R_wave_limits.V5.R_wave_end);
V6_start = round(R_wave_limits.V6.R_wave_start);
V6_end = round(R_wave_limits.V6.R_wave_end);

X_R_wave = [ecg.I(I_start:I_end) ecg.II(II_start:II_end) ...
    ecg.V1(V1_start:V1_end) ecg.V2(V2_start:V2_end) ...
    ecg.V3(V3_start:V3_end) ecg.V4(V4_start:V4_end) ...
    ecg.V5(V5_start:V5_end) ecg.V6(V6_start:V6_end)]';

% Cut out T wave

X = [ecg.I(fidpts(3):fidpts(4)) ecg.II(fidpts(3):fidpts(4)) ...
    ecg.V1(fidpts(3):fidpts(4)) ecg.V2(fidpts(3):fidpts(4)) ...
    ecg.V3(fidpts(3):fidpts(4)) ecg.V4(fidpts(3):fidpts(4)) ...
    ecg.V5(fidpts(3):fidpts(4)) ecg.V6(fidpts(3):fidpts(4))]';

% DONT Subtract out the centroid to avoid offsets since already zerod data

% SVD - Only care about left singular vectors U and singular values, but
% easier to extract singular values using second line of code than to get
% the entire matrix S with all the zeros etc.
[U,~,~]=svd(X);
s=svd(X);
s_R_wave=svd(X_R_wave);


%%% TMD

% See Concept of T-Wave Morphology Dispersion
% https://ieeexplore.ieee.org/document/825905

% Take first 2 singular vectors
Ut = U(:,[1 2]);

% Make 2x2 matrix of singular values on diagonals
S = [s(1) 0 ; 0 s(2)];

% Create W
W = (Ut * S)';

% Drop V1 (3rd col of W)
W(:,3) = [];

% Calculate angles between pairs of columns of W

% Preallocate A to store permutations of angles.  Since have 7 leads after
% excluding V1 it will be 7x7
A = zeros(7,7);

% Calculate angle permutations (set z = 0 to use cross function)
for i = 1:7
    for j = 1:7
        wi = [W(:,i) ; 0];
        wj = [W(:,j) ; 0];
        A(i,j) = atan2d(norm(cross(wi,wj)),dot(wi,wj));
    end
end

% Want to exclude if i = j
% Take values above diagnoal; there are 21 values
A2 = triu(A,1);

% Average the 21 values in A2 to get TMD
TMD = sum(A2,"all")/21;


%%% TWR

% See Analysis of T wave morphology parameters with signal averaging during 
% ischemia induced by percutaneous transluminal coronary angioplasty
% https://ieeexplore.ieee.org/document/5445289
%
% See Role of Dipolar and Nondipolar Components of the T Wave in Determining 
% the T Wave Residuum in an Isolated Rabbit Heart Model
% https://onlinelibrary.wiley.com/doi/10.1046/j.1540-8167.2004.03466.x
%
% TWR is sum of 4th through 8th eigenvalues 
% (although some other reports say that they use the squares of the eigenvalues?)
% The singular values are the sqrts of the eigenvalues
% Therefore the eigenvalues are the singular values squared
%
% The sum of the 4th to 8th eigenvalues expressed in mV2 defines the 
% TWR and quantifies the nondipolar components of the T wave.

% Generate eigenvalues as 'eigmat' by squaring singular values
eigmat = s_R_wave.^2

% TWR is sum of 4th through nth (8th in this case) eigenalues
RWR_abs = sum(eigmat(4:8))
fprintf(RWR_abs)
% Relative TWR is the percentage of the whole
RWR_rel = 100 * TWR_abs / sum(eigmat)
fprintf(RWR_rel)
end



function R_wave_limits = find_R_wave_start_end_and_average(ecg, fidpts)
    % Function to find R_wave_start and R_wave_end for each lead (I, II, V1-V6),
    % and then average these points across all leads.
    % Inputs:
    %   ecg - structure containing lead signals (e.g., ecg.I, ecg.II, ecg.V1, etc.)
    %   fidpts - array with fiducial points; fidpts(1) is QRS start, fidpts(2) is R peak, fidpts(3) is S end
    % Output:
    %   R_wave_limits - structure with individual and averaged R_wave_start and R_wave_end

    % Define leads to analyze
    lead_names = {'I', 'II', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
    
    % Get fiducial points
    QRS_start = fidpts(1);
    R_peak = fidpts(2);
    S_end = fidpts(3);
    
    % Initialize arrays to store R wave start and end values for averaging
    R_wave_start_values = NaN(1, length(lead_names)); % NaN array for start averaging
    R_wave_end_values = NaN(1, length(lead_names));   % NaN array for end averaging

    % Loop over each lead to calculate R_wave_start and R_wave_end
    for i = 1:length(lead_names)
        lead = lead_names{i};
        
        % Initialize R_wave_start and R_wave_end as empty
        R_wave_start = [];
        R_wave_end = [];

        % Extract the data for the current lead
        lead_data = ecg.(lead);

        %% Find R_wave_start by searching backward from R_peak
        for j = R_peak:-1:QRS_start+1
            % Check for sign change between consecutive points
            if sign(lead_data(j)) ~= sign(lead_data(j - 1))
                R_wave_start = j; % Set R_wave_start to the first sign change
                break;
            end
        end

        %% Find R_wave_end by searching forward from R_peak
        for j = R_peak:S_end-1
            % Check for sign change between consecutive points
            if sign(lead_data(j)) ~= sign(lead_data(j + 1))
                R_wave_end = j + 1; % Set R_wave_end to the point after the sign change
                break;
            end
        end

        % Store R_wave_start and R_wave_end in the output structure and arrays for averaging
        R_wave_limits.(lead).R_wave_start = R_wave_start;
        R_wave_limits.(lead).R_wave_end = R_wave_end;

        if ~isempty(R_wave_start)
            R_wave_start_values(i) = R_wave_start; % Save start point if found
        end

        if ~isempty(R_wave_end)
            R_wave_end_values(i) = R_wave_end; % Save end point if found
        end
    end

    %% Calculate Averaged Start and End Points Across All Leads
    % Take the mean of all detected start and end points, omitting NaN values
    avg_start = round(mean(R_wave_start_values, 'omitnan'));
    avg_end = round(mean(R_wave_end_values, 'omitnan'));

    % Store the averaged results in the output structure
    R_wave_limits.avg_start = avg_start;
    R_wave_limits.avg_end = avg_end;

    % Display the results
    disp('R_wave_start and R_wave_end for each lead:');
    disp(R_wave_limits);
    fprintf('Averaged R_wave_start: %d, Averaged R_wave_end: %d\n', avg_start, avg_end);
end
