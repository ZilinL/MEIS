%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paper: Manifold Embedded Instance Selection to Suppress Negative Transfer in Motor Imagery-based Brain-Computer Interface
%% Author: Code by Zilin Liang
%% Date: 2023/Oct./9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear; close all; warning off;

% Load datasets:
root='data\MI1\';
listing=dir([root '*.mat']);
addpath('lib');

fnum=length(listing);
Ca=nan(59,59,200*fnum);
Xr=nan(59,300,200*7);
Xa=nan(59,300,200*7);
Y=nan(200*fnum,1);
ref={'riemann','logeuclid','euclid'};
for f=1:fnum
    load([root listing(f).name])
    idf=(f-1)*200+1:f*200;
    Y(idf) = y; Xr(:,:,idf) = x;
    Ca(:,:,idf) = centroid_align(x,ref{1});
    [~,Xa(:,:,idf)] = centroid_align(x,ref{3});
end

%%%%%%%%%%%%%%%%%according to Riemann geodesic distance%%%%%%%%%%%%%%%%%%%%
% By calculating the Riemann geodesic distance, the subject number of the source domain closest to each target domain is obtained.
HSS_sub = {6,6,7,2,2,1,1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BCA=zeros(fnum,1);
BCA_rw=zeros(fnum,1);
for n=1:fnum
    ids=[];
    disp(n)
    % Single target data & single source data
    idt=(n-1)*200+1:n*200;
    HSS_s_sub = HSS_sub{n};
    num_HSS_s_sub = length(HSS_sub{n});
    for i=1:num_HSS_s_sub
        n_sub = HSS_sub{n}(i);
        ids = [ids (n_sub-1)*200+1:n_sub*200];
    end
    Yt=Y(idt); Ys=Y(ids);
    idsP=Yt==1; idsN=Yt==0;
    Ct=Ca(:,:,idt);  Cs=Ca(:,:,ids);

    % Logarithmic mapping on aligned covariance matrices
    Xs=logmap(Cs,'MI'); % dimension: 253*1152 (features*samples)
    Xt=logmap(Ct,'MI');

    % Dimensionality reduction by one-way ANOVA based on F-values
    [idx, Fs]=Fvalue(Xs', Ys, length(Ys));
    Xs=Fs'; Xt=Xt(idx,:);

    %% MEIS
    options.d = 8;             % subspace bases
    options.T = 5;              % iterations, default=5
    options.alpha= 0.01;        % the parameter for source discriminability, default=0.01, max=0.0001
    options.beta = 1;         % the parameter for target locality, default=0.1
    options.rho = 800;           % the parameter for subspace discrepancy
    options.clf = 'slda';        % the string for base classifier, 'slda' or 'svm'
    Cls = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%All samples%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Zs, Zt, rw] = MEIS(Xs, Xt, Ys, Cls, options);
    Ypre = slda(Zt,Zs,Ys);
    BCA(n)=.5*(mean(Ypre(idsP)==1)+mean(Ypre(idsN)==0));

    %%%%%%%%%%%%%%Remove 20% of negative migration samples%%%%%%%%%%%%%%%%%
    [B,I] = sort(rw,'descend');
    nI = floor(length(I)*0.1);
    Zs_rw = Zs(:,I(nI+1:end-nI));
    Ys_rw = Ys(I(nI+1:end-nI));
    Ypre = slda(Zt,Zs_rw,Ys_rw);
    BCA_rw(n)=.5*(mean(Ypre(idsP)==1)+mean(Ypre(idsN)==0));

end
disp(['BCA:',num2str(mean(BCA)*100)]);
disp(['BCA_rw:',num2str(mean(BCA_rw)*100)]);

