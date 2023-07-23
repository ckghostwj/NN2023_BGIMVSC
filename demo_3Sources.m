% This is the demo code for the paper:
% Lilei Sun, Jie Wen, et al., Balance guided incomplete multi-view spectral clustering, Neural Networks, 2023.
% For any problems, please contact Lilei Sun via sunlileisun@163.com  or Jie Wen via jiwen_pr@126.com.
% If you use the code, please cite the above reference.
% For different version of Matlab, the result may be different.

clear all;
clc

Dataname = '3sources3vbigRnSp';
Datafold = '3sources3vbigRnSp_missing_rate_0.5.mat';

load(Dataname);
load(Datafold);
[numFold,numInst] = size(folds);
numClust = length(unique(truth));
numInst  = length(truth);

para_r = 7;
lambda = 0.1;
gamma = 10;
miu = 0.00001;
rho = 1.05;
k = 5;
rand('seed',5867);
f = 2;

load(Dataname);
ind_folds = folds{f};
truthF = truth;
if size(X{1},2)~=length(truth)
    for iv = 1:length(X)
        X{iv} = X{iv}';
    end
end
linshi_AAW = 0;
linshi_WWW= 0;
S_ini = 0;
for iv = 1:length(X)
    X1 = X{iv};
    X1 = NormalizeFea(X1,0);
    aa = sum(ind_folds,2);
    ind_1 = find(ind_folds(:,iv) == 1);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(:,ind_0) = [];
    Y{iv} = X1;
    linshi_W = ones(numInst,numInst);
    linshi_W(:,ind_0) = 0;
    linshi_W(ind_0,:) = 0;
    Wt{iv} = linshi_W;
    options = [];
    options.NeighborMode = 'KNN';
    options.k = k;
    options.WeightMode = 'Binary';
    So = constructW(X1',options);
    
    G = diag(ind_folds(:,iv));
    G(:,ind_0) = [];
    Sor{iv} = G*So*G';
    linshi_AAW = linshi_AAW+Sor{iv}.*linshi_W;
    linshi_WWW = linshi_WWW+linshi_W;
    S_ini = S_ini+Sor{iv};
end
clear X X1 W1 ind_0 So
X = Y;
clear Y
max_iter = 100;
[P,S,Q,obj] = V11_Function(S_ini,Sor,Wt,numClust,lambda,gamma,miu,rho,para_r,max_iter);

[~,pre_labels] = max(P',[],1);
result_cluster = ClusteringMeasure(truthF, pre_labels)*100;
acc(f) = result_cluster(1);
nmi(f) = result_cluster(2);
pur(f) = result_cluster(3);
if size(pre_labels,1)~=size(truthF,1)
    pre_labels = pre_labels';
end

pre_labels_kmeans    = kmeans(P,numClust,'emptyaction','singleton','replicates',20,'display','off');
result_kmeans = ClusteringMeasure(truthF,pre_labels_kmeans)*100
