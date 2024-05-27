% Alternating minimization for Dantzig Selector Model

%%==== Input variables ========
%         X: Design matrix
%         y: observation vertor
%     delta: tuning parameter
%      beta: the vector of regression coefficients
%  fixp.eps: stopping precision
%  fixp.MAX: the maximum iteration
%         D: diagonal matrix whose diagonal entries are the norm of the
%            columns of X;

% Corresponding to Hongjin He:
% Email to : hehjmath@hdu.edu.cn

function out = PLAM(X,D,y,delta,para,fixp)
[~,p] = size(X);     b = (y'*X)';    ddel = delta.*D;  beta = zeros(p,1);
z = X'*(X*beta) - b;
gamma = para.gamma;                      % the penalty parameter, gamma >0;
mu = para.mu;                            % the proximal parameter;
gm = gamma*mu;                           % gamma * mu

tic;     id = 0;

for iter = 1 : fixp.MAX
    XbX = X'*(X*beta);
    temp =  XbX - b ;
    znew = max( min( temp , ddel), -ddel);              % Projection onto infty-norm box    
    beta_new = shrink( beta - ( X'*(X*(temp-znew)) )./mu, 1/gm); % shrinkage
    
    e1 = norm(beta - beta_new, inf) / max( norm(beta_new), 1 ) ;
    e2 = norm(z - znew,inf) / max( norm(znew),1 );
    out.error(iter) = max( e1, e2 );
    if fixp.detail == 1
        fprintf('Iter = %d && Err1 = %2.5e && Err2 = %2.5e && Error = %2.5e \n ',iter,e1,e2,out.error(end))
    end
    
    beta = beta_new;    z = znew;
    
    if iter == 1 || mod(iter,5) == 0 || out.error(end) <= fixp.eps
        id = id + 1;
        out.objs(id) = norm(beta,1);   out.iters(id) = iter;
    end
    
    
    if out.error(end) <= fixp.eps || iter >= fixp.MAX
        out.iter = iter;  out.solution = beta;  out.obj = out.objs(end);
        break;
    end
    
end
out.time = toc;
end

% shrinkage operator; sol = argmin \mu*||x||_1 + 0.5*||x -y||^2
function sol = shrink(y,mu)
sol = sign(y).*max( 0, abs(y) - mu );
end
