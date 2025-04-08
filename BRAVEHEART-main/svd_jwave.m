function [JWR_abs, JWR_rel] = svd_jwave(ecg, J_wave_limits)
start_pt = J_wave_limits.start;
% end_pt = J_wave_limits.end;
end_pt = start_pt+40;
% if end_pt < start_pt+8
%     end_pt = start_pt+8;
% end
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
JWR_abs = sum(eigmat(4:8));

% Relative TWR is the percentage of the whole
JWR_rel = 100 * JWR_abs / sum(eigmat);
