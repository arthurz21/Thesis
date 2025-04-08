%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TMD_i = svd_twave_2(ecg, pt1, pt2)


% Cut out T wave
X = 1000*[ecg.I(pt1:pt2) ecg.II(pt1:pt2) ...
    ecg.V1(pt1:pt2) ecg.V2(pt1:pt2) ...
    ecg.V3(pt1:pt2) ecg.V4(pt1:pt2) ...
    ecg.V5(pt1:pt2) ecg.V6(pt1:pt2)]';

% DONT Subtract out the centroid to avoid offsets since already zerod data

% SVD - Only care about left singular vectors U and singular values, but
% easier to extract singular values using second line of code than to get
% the entire matrix S with all the zeros etc.
[U,~,~]=svd(X);
s=svd(X);


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
TMD_i = sum(A2,"all")/21;
