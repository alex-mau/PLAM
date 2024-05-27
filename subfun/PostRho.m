function  [xp,rho1,rho2] = PostRho(x,X,y,betas,sigma,peps)
% ===== Output variables:
%  xp : The postprocessing solution after post processing
% rho1: The original rho--- The preprocessing errors
% rho2: The new rho--- The postprocessing error;

% ===== input variables:
%     x: the preprocessing solution obtained via the method;
%     X: the design matrix
%     y: The observation
% betas: the optimal solution
% sigma: the standard deviation of noise
%  peps: the precision of postprocessing;
% if class(x)=='gpuArray'
%     x=gather(x);
%     X=gather(X);
%     y=gather(y);
%     betas=gather(betas);
% end
xp = PostPro(x,X,y,peps);               % Postprocessing
rho1 = norm(x - betas)^2 / sum( min(betas.^2,sigma^2) );
rho2 = norm(xp - betas)^2 / sum( min(betas.^2,sigma^2) );
