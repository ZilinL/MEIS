  
function [srk] = sample_reweight_kmm(X, Z, dis, Theta)
% Author: Zilin Liang
% Paper: Manifold Embedded Instance Selection to Suppress Negative Transfer in Motor Imagery-based Brain-Computer Interface
% Date: 2023/Oct./9
% Use Kernel Mean Matching to estimate weights for importance weighting.
% Code modification from the article "Correcting Sample Selection Bias by unlabeled data"

% IPUTS:
% X: Source domain data.
% Z: Target domain data.
% dis: Distance metric method, can be 'Euclidean', 'AIRM', 'Stein', 'Jeffrey', or 'Log-Euclidean'.
% Theta: Parameter controlling the kernel width.

% OUTPUTS:
% srk: Sample weights, size 1 x NX, used to adjust the weights of source domain data.

% Parse input
theta =Theta;
kernel = 'rbf';
distance = dis;

switch kernel
    case 'rbf'
        switch distance
            case 'Euclidean'
                % Shapes
                [NX,~] = size(X);
                [NZ,~] = size(Z);
                epsilon = (sqrt(NX)-1)/sqrt(NX);
                % Calculate Euclidean distances
                K = pdist2(X, X);
                k = pdist2(X, Z);
            case 'AIRM'
                % Shapes
                NX = size(X,3);
                NZ= size(Z,3);
                epsilon = 1/sqrt(NX);
                K = compute_distance(X,X,'A');
                k = compute_distance(X,Z,'A');
            case 'Stein'
                % distance while using Stein
                K = compute_distance(X,X,'S');
                k = compute_distance(X,Z,'S');
            case 'Jeffrey'
                % distance while using Jeffrey
                K = compute_distance(X,X,'J');
                k = compute_distance(X,Z,'J');
            case 'Log-Euclidean'
                % distance while using LogED
                K = compute_distance(X,X,'L');
                k = compute_distance(X,Z,'L');
        end
        
        
        % Cleanup
        I = find(K<0); K(I) = zeros(size(I));
        J = find(K<0); K(J) = zeros(size(J));
        
        % Radial basis function
        K = exp(-K/(theta));
        k = exp(-k/(theta));
        k = 2*NX./NZ*sum(k,2);
        
    case 'diste'
        % Calculate Euclidean distances
        K = pdist2(X, X);
        k = pdist2(X, Z);
        if theta ~= 2
            K = sqrt(K).^theta;
            k = sqrt(k).^theta;
        end
        k = 2*NX./NZ*sum(k,2);
        
end

% Solve quadratic program
options.Display = 'off';
 srk = quadprog(K,-k,[ones(1,NX); -ones(1,NX)],[NX.*epsilon+NX NX.*epsilon-NX],[],[],zeros(NX,1),3*ones(NX,1), [], options)';
end