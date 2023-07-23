function [P,S,Q,obj] = V11_Function(S_ini,Sor,Wt,numClust,lambda,gamma,miu,rho,para_r,max_iter)

% ------- ³õÊ¼»¯ ------- %
S = S_ini;
S = (S+S')*0.5;
LSv = diag(sum(S,1))-S;
[F, ~, ~] = eig1(LSv, numClust, 0);
clear LSv

P = F;
P = max(P,0);
Q = P;
C = zeros(size(P));
alpha = ones(length(Sor),1)/length(Sor);
alpha_r = alpha.^para_r;

for iter = 1:max_iter
    
    
    % --------- S --------- %
    linshi_AAW = 0;
    linshi_WWW= 0;
    for iv = 1:length(Sor)
        linshi_AAW = linshi_AAW+alpha_r(iv)*Sor{iv}.*Wt{iv};
        linshi_WWW = linshi_WWW+alpha_r(iv)*Wt{iv}; 
    end
    G = P*F';
    H = (lambda*G+linshi_AAW)./(lambda+linshi_WWW);
    H = H-diag(diag(H));
    S = zeros(size(H));
    
    for is = 1:size(H,1)
        ic = [1:size(H,1)];
        ic(is) = [];
        tmp_s = EProjSimplex_new(H(ic,is));
        S(ic,is) = tmp_s;
    end
    % --------- F ------- %
    W = S'*P;
    W(isnan(W)) = 0;
    W(isinf(W)) = 0;
    [U,~,V] = svd(W,'econ');
    U(isnan(U)) = 0;
    U(isinf(U)) = 0;
    V(isnan(V)) = 0;
    V(isinf(V)) = 0;
    F = U*V';
    % --------- P --------- %
    linshi_P = (2*lambda*S*F+miu*Q-C)/(2*lambda+miu);
    P = zeros(size(linshi_P));
    
    for is = 1:size(linshi_P,1)
        P(is,:) = EProjSimplex_new(linshi_P(is,:));
    end
    % --------- Q --------- %
%     Q = inv(gamma * ones(size(S, 1), size(S, 1)) + miu * eye(size(S, 1))) * (miu * P + C);
    
    Nsamp = size(S,1);
    Q = (1/(miu*(miu+Nsamp*gamma)))*((miu+Nsamp*gamma)*eye(Nsamp,Nsamp)-gamma*ones(Nsamp,Nsamp))* (miu * P + C);
    % --------- alpha -------- %
    for iv = 1:length(Sor)
        Rec_error(iv) = norm((S-Sor{iv}).*Wt{iv},'fro')^2;
    end
    aH = bsxfun(@power,Rec_error, 1/(1-para_r));     % h = h.^(1/(1-r));
    alpha = bsxfun(@rdivide,aH,sum(aH)); % alpha = H./sum(H);
    alpha_r = alpha.^para_r;
    % --------- C --------- %
    
    C = C + miu * (P - Q);
    miu = min(miu*rho,1e8);
    % ----- obj ------ %
    obj_iter = alpha_r*Rec_error'+lambda*norm(S-P*F')+0.5*gamma*trace(P'*ones(size(P,1),size(P,1))*P);
    obj(iter) =  obj_iter;
    if iter > 2 && abs(obj(iter)-obj(iter-1))<1e-5
        iter;
        break;
    end
    leq = P-Q;
    if max(abs(leq)) < 1e-3
        %fprintf('leq');
        iter;
        break;
    end
end
end