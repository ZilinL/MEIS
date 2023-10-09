function [Zs, Zt, srk] = MEIS(Xs, Xt, Ys, Yt0, options)
% Author: Zilin Liang
% Paper: Manifold Embedded Instance Selection to Suppress Negative Transfer in Motor Imagery-based Brain-Computer Interface
% Date: 2023/Oct./9
%
% Input: Xs: Source domain features, in D*N format, where D is feature dimension and N is the number of samples.
%        Xt: Target domain features, in D*M format, where D is feature dimension and M is the number of samples.
%        Ys: Source domain labels, 1*N.
%        Yt0: Target pseudo-labels.
% Options - Struct containing parameters:
%        d: Subspace dimension.
%        T: Number of iterations, default=5.
%        beta: the parameter for U,
%        alpha: the parameter for P,
%        rho: the parameter for Q,
%        clf: the string for base classifier, 'slda' or 'svm'.
% Output: 
% Zs: Embeddings for source domain features, in d*N format.
% Zt: Embeddings for target domain features, in d*M format.
% srk: Weights for source domain samples, an array of size 1*N, used for sample reweighting.

% Set options
d = options.d; T = options.T;
alpha = options.alpha;
beta = options.beta;
rho = options.rho; clf = options.clf;

% Get variable sizes
[ms,ns] = size(Xs); [mt,nt] = size(Xt);
class = unique(Ys); C = length(class);


% Initialize srk
Yt0 = slda(Xt,Xs,Ys);
srk = ones(ns,1);
MaxIter = 20;
A = eye(ms); B = eye(ms);

for iMaxIter = 1:MaxIter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% weight calculation %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = 80;
    [srk] = sample_reweight_kmm(Xs'*A, Xt'*B, 'Euclidean', theta);
    srk = srk';
    con(iMaxIter) = norm(srk, 'fro')^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% matrix W calculation %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % U: Minimize Source Domain Intra-Class Variance and Maximize Inter-Class Variance
    [Ssw, Ssb] = ScatterMatrix(Xs, Ys);
    U = zeros(2*ms,2*ms); U(1:ms,1:ms) = Ssw;
    U0 = zeros(2*ms,2*ms); U0(1:ms,1:ms) = Ssb;

    % S: Target Domain Variance Maximization
    Ht = eye(nt)-1/(nt)*ones(nt,nt);
    P = [zeros(ms),zeros(mt); zeros(ms),Xt*Ht*Xt'];

    % K: Subspace Discrepancy Minimization |B-A|_F+|B|_F
    Q = [eye(ms),-eye(mt);-eye(ms),2*eye(mt)];

    % V£ºMinimization of the local projection map
    manifold.k = 10; % default set to 10
    manifold.NeighborMode = 'KNN';
    manifold.WeightMode = 'HeatKernel';
    S = lapgraph(Xt',manifold);
    D = full(diag(sum(S,2)));
    L = D-S;
    L = [zeros(ms),zeros(mt); zeros(ms),Xt*L*Xt'];

    for t = 1:T

        % Z: Maximum Mean Discrepancy (MMD) Distance Minimization
        Xs_rw = Xs*diag(srk);
        e0 = [1/ns*srk;-1/nt*ones(nt,1)];
        M0 = e0*e0';
        Mc = 0;
        for c = reshape(unique(Ys),1,length(unique(Ys)))
            ec = zeros(ns + nt,1);
            ec(Ys == c) = 1 / length(find(Ys == c));
            ec(ns + find(Yt0 == c)) = -1 / length(find(Yt0 == c));
            ec(isinf(ec)) = 0;
            Mc = Mc + ec * ec';
        end
        M = M0 + Mc;
        M = M / norm(M,'fro');
        X = [Xs_rw,zeros(size(Xt));zeros(size(Xs_rw)),Xt];
        Z = X*M*X';

        % generalized eigen-decomposition
        Emin = alpha*U +  beta*L + rho*Q + Z;
        Emax = alpha*U0 + P; % Matrix V
        [W,~] = eigs(Emin+10^(-3)*eye(ms+mt), Emax, d, 'SM'); % SM: smallestabs

        % Smallest magnitudes
        A = W(1:ms, :);
        B = W(ms+1:end, :);

        % Embeddings
        Zs = A'*Xs;
        Zt = B'*Xt;

        if t>1
            if strcmp(clf,'slda')
                Yt0 = slda(Zt,Zs,Ys);
            elseif strcmp(clf,'svm')
                w=ones(size(Ys)); w(Ys==1)=sum(Ys==0)/sum(Ys==1);
                model = libsvmtrain(w,Ys,Zs','-h 0 -t 0 -c 0.125 -q');
                Yt0 = libsvmpredict(ones(size(Zt,2),1),Zt',model);
            end
        end
    end
    con_A(iMaxIter) = norm(A, 'fro')^2;
    con_B(iMaxIter) = norm(B, 'fro')^2;


    %%%%%%%%%%%%%%iterating the objective function%%%%%%%%%%%%%%%
    con_ObjFunc(iMaxIter) = norm(Z, 'fro')^2+alpha*trace(A'*Ssw*A)+beta*trace(B'*Xt*(D-S)*Xt'*B)+rho*trace(norm(A-B, 'fro')^2+norm(B, 'fro')^2);

end

end

function y_onehot=onehot(y,class)
% Encode label to onehot form
% Input:
% y: label vector, N*1
% Output:
% y_onehot: onehot label matrix, N*C

nc=length(class);
num=length(y);
y_onehot=zeros(num, nc);
for i=1:num
    idx=nc-find(class==y(i))+1;
    y_onehot(i, idx)=1;
end
end

function [Sw,Sb] = ScatterMatrix(X, Y)
[m,n] = size(X);
Sw = zeros(m);Sb = zeros(m);
class = unique(Y);C = length(class);
meanTotal = mean(X,2);
for i=1:C
    Xi = X(:,Y==class(i));
    meanClass = mean(Xi,2);
    Hi = eye(size(Xi,2))-1/(size(Xi,2))*ones(size(Xi,2),size(Xi,2));
    Sw = Sw + Xi*Hi*Xi';
    Sb = Sb + size(Xi,2)*(meanClass-meanTotal)*(meanClass-meanTotal)';
end

end
