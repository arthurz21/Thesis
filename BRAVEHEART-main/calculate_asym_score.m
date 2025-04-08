function [Tasym] = calculate_asym_score(vcg, T_start, T_end)
% Example usage:
% Assuming you have your VCG data and T wave boundaries:
x = vcg.X; % X lead of VCG
y = vcg.Y; % Y lead of VCG
z = vcg.Z; % Z lead of VCG
pca1 = calculatePCA1(x, y, z, T_start, T_end);
T_start = 1;
T_end = length(pca1);
T_wave_seg = pca1(T_start:T_end);
[max_val, local_peak_idx] = max(T_wave_seg);
t_peak_idx = T_start + local_peak_idx - 1;

% Calculate first derivative of T wave
derivative = diff(T_wave_seg);
derivative = [derivative; derivative(end)];

% Split into ascending and descending segments
ascending = derivative(1:local_peak_idx);
descending = derivative(local_peak_idx+1:end);

% Check for empty segments first
if isempty(ascending) || isempty(descending)
    warning('Empty segment detected in T wave - peak may be at edge of wave');
    Tasym = NaN;
    return;
end

% Check if segments are too short (need at least 2 points for interpolation)
if length(ascending) < 2 || length(descending) < 2
    warning('Segments too short for analysis - peak may be at edge of wave');
    Tasym = NaN;
    return;
end

% Normalize each segment by its maximum absolute amplitude
asc_max = max(abs(ascending));
desc_max = max(abs(descending));

if asc_max == 0 || desc_max == 0 || isnan(asc_max) || isnan(desc_max)
    warning('Zero amplitude or NaN detected in T wave segment');
    Tasym = NaN;
    return;
end

ascending_norm = ascending / asc_max;
descending_norm = descending / desc_max;

% Flip descending segment (both x and y axis)
descending_flipped = -flipud(descending_norm);

% Interpolate to match lengths if necessary
if length(ascending_norm) ~= length(descending_flipped)
    % Determine target length - use the longer of the two segments
    target_length = max(length(ascending_norm), length(descending_flipped));
    
    % Use proper indexing for interpolation
    orig_asc_indices = 1:length(ascending_norm);
    orig_desc_indices = 1:length(descending_flipped);
    
    % Target indices for interpolation
    target_indices = linspace(1, length(ascending_norm), target_length);
    ascending_norm = interp1(orig_asc_indices, ascending_norm, target_indices, 'linear', 'extrap');
    
    target_indices = linspace(1, length(descending_flipped), target_length);
    descending_flipped = interp1(orig_desc_indices, descending_flipped, target_indices, 'linear', 'extrap');
end

% Calculate residuals
r = ascending_norm - descending_flipped;
n = length(r);
Tasym = sum(r.^2)/(n-1);
end

function pca1 = calculatePCA1(x_lead, y_lead, z_lead, t_start, t_end)
% Calculate PCA1 from VCG XYZ leads
% Inputs:
% x_lead, y_lead, z_lead - Full VCG XYZ leads
% t_start - T wave start index
% t_end - T wave end index

% Extract ST-T segment from each lead
x_t = x_lead(t_start:t_end);
y_t = y_lead(t_start:t_end);
z_t = z_lead(t_start:t_end);

% Create data matrix where each column is a lead
X = [x_t, y_t, z_t];

% Center the data (subtract mean)
X = X - mean(X);

% Perform PCA using SVD
[U, S, V] = svd(X, 'econ');

% First principal component scores
pca1 = X * V(:,1);
end