% Partially Linearized ADMM for Dantzig selector;
% Refer to: Regularized Split Feasibility Problem
% Authors: Hongjin He, Bingsheng HE, and Hongkun XU

%%==== Input variables ========
%         X: Design matrix: Unit coloumn
%         y: observation vertor
%     delta: tuning parameter
%      beta: the vector of regression coefficients
%  fixp.eps: stopping precision
%  fixp.MAX: the maximum iteration
%      diag: diagonal matrix whose diagonal entries are the norm of the
%            columns of X;
%     Mtype: The type of your method, Original PLADM or PLADM with
%            an extra correction step

% Corresponding to Hongjin He:
% Email to : hehjmath@hdu.edu.cn

function out = PLADM(X,diag,y,delta,Mtype,para,fixp)

[~,p] = size(X);     b = X'*y;   ddel = delta.*diag;
beta = zeros(p,1);   lam = beta;         %initial points;

gamma = para.gamma;                      % the penalty parameter, gamma >0;
mu = para.mu;                            % the proximal parameter;
gm = gamma/mu;

if strcmp(Mtype,'correction')
    tau = para.tau;               % the relaxtion factor for acceleration
end

tic;   XbX = X'*(X*beta);  out.time = toc;  id = 0;

for iter = 1 : fixp.MAX
    tic;
    temp =  XbX - b - lam./gamma ;
    znew = max( min( temp , ddel), -ddel);              % Projection onto Q
    
    beta_temp = shrink( beta - gm.*( X'*(X*(temp-znew)) ), 1/mu); % shrinkage
    
    XbXt =  X'*(X*beta_temp);
    lam_temp = lam - gamma*( XbXt- b - znew );
    
    if strcmp(Mtype,'correction')
        ldl = lam_temp - lam;  bdb = beta_temp - beta;
        domi_a = (bdb'*bdb)*mu + (ldl'*ldl)/gamma;
        domi_b = ldl'* ( XbXt - XbX );
        alpha = 1 + domi_b/domi_a ;
        beta_new = beta + tau*alpha.*bdb;
        lam_new = lam + tau*alpha.*ldl;
    else
        beta_new = beta_temp;
        lam_new = lam_temp;
    end
    XbX = X'*(X*beta_new);
    
    if strcmp(fixp.rule,'SRI')
        e1 = norm(beta - beta_new, inf) / max( norm(beta_new), 1 ) ;
        e2 = norm(lam - lam_new,inf) / max( norm([beta_new;lam_new]),1 );
        out.error = max( e1, e2 );
    elseif strcmp(fixp.rule,'SRII')
        e1 = norm(beta - beta_new, inf) / max( norm(beta_new), 10 ) ;
        e2 = norm(lam - lam_new,inf) / (10*norm([beta_new;lam_new]));
        out.error = max( e1, e2 );
    elseif strcmp(fixp.rule,'SRIII')  
        e1 = max( abs(beta - beta_new) ) / max( norm(beta_new), 1 ) ;
        e2 = (norm( ( XbX - b)./diag,inf) - delta)/max(norm(beta_new),1);
        out.error = max(e1,e2);
    else
        display('Please refer to PLADM.m !!')
    end
    
    beta = beta_new;  lam = lam_new;       t = toc;
    
    out.time = out.time + t;
    
    if iter == 1 || mod(iter,5) == 0 || out.error <= fixp.eps
        id = id + 1;
        out.objs(id) = norm(beta,1);
        out.iters(id) = iter;
    end
    if out.error <= fixp.eps || iter >= fixp.MAX
        out.iter = iter;  out.beta = beta;  out.obj = out.objs(end);
        break;
    end
    if fixp.detail == 1
        fprintf('Iter = %d && Err1 = %2.5e && Err2 = %2.5e  \n ',iter,e1,e2)
    end
    
end

end

% shrinkage operator; sol = argmin \mu*||x||_1 + 0.5*||x -y||^2
function sol = shrink(y,mu)
sol = sign(y).*max( 0, abs(y) - mu );
end
