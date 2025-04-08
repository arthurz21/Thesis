function [STR_abs, STR_rel] = svd_stwave(ecg, ST_segment_limits)

start_pt = ST_segment_limits.ST_start;
end_pt = ST_segment_limits.ST_end;
if isnan(end_pt)
    end_pt = start_pt +8;
end
if end_pt < start_pt +8
    end_pt = start_pt +8;
end 
% Cut out R wave

X = [ecg.I(start_pt:end_pt) ecg.II(start_pt:end_pt) ...
    ecg.V1(start_pt:end_pt) ecg.V2(start_pt:end_pt) ...
    ecg.V3(start_pt:end_pt) ecg.V4(start_pt:end_pt) ...
    ecg.V5(start_pt:end_pt) ecg.V6(start_pt:end_pt)]';

% DONT Subtract out the centroid to avoid offsets since already zerod data

% SVD - Only care about left singular vectors U and singular values, but
% easier to extract singular values using second line of code than to get
% the entire matrix S with all the zeros etc.
s=svd(X);
eigmat = s.^2;

% TWR is sum of 4th through nth (8th in this case) eigenalues
STR_abs = sum(eigmat(4:8));

% Relative TWR is the percentage of the whole
STR_rel = 100 * STR_abs / sum(eigmat);
end