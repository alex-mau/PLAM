% Customized Proximal Point Algorithm for Dantzig selector;
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
%     Mtype: The type of your method, Dual-Primal or Primal-Dual 

% Corresponding to Hongjin He:
% Email to : hehjmath@hdu.edu.cn

function out = CPPA(X,diag,y,delta,Mtype,para,fixp)

[~,p] = size(X);        b = X'*y;   ddel = delta.*diag;
 beta = zeros(p,1);   lam = beta;             %initial points;

r = para.r;                % the parameter associated with beta
s = para.s;                % the parameter associated with Lagrangian lam;
tau = para.tau;            % the relaxtion factor tau \in (0,2);

out.time = 0;  id = 0;

for iter = 1 : fixp.MAX
    tic;
    if strcmp(Mtype,'DPM')
        temp = X'*(X*beta) - b - s.*lam ;
        lam_temp = ( max( min( temp , ddel), -ddel) - temp)./s ;             % Dual step
        beta_temp = shrink( beta + ( X'*(X*(2.*lam_temp - lam)) )./r , 1/r); % Primal step
    elseif strcmp(Mtype,'PDM')
        beta_temp = shrink( beta + (X'*(X*lam))./r, 1/r);    % Primal step
        temp = X'*(X*(2.*beta_temp - beta)) - b - s.*lam;
        lam_temp = ( max( min( temp , ddel), -ddel) - temp)./s ;
    else
        display('Please refer to CPPA.m!!');
    end
    beta_new = beta - tau*( beta - beta_temp);  % Update the new variable
    lam_new  = lam  - tau*( lam - lam_temp);    % Update the multiplier
   
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
        e2 = (norm( ( XbX- b)./diag,inf) - delta)/max(norm(beta_new),1);
        out.error = max(e1,e2);
    else
        display('Please refer to CPPA.m !!')
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
    
end

end

% shrinkage operator; sol = argmin \mu*||x||_1 + 0.5*||x -y||^2
function sol = shrink(y,mu)
sol = sign(y).*max( 0, abs(y) - mu );
end
