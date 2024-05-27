%  The nonmontone gradient method appoximately solves
%
%   min ||x||_1 + (mu/2)*||A'Ax - b||^2.
%
function [x,iter_in] = spg2(A,Atype,b,mu,x,eps,maxiter)

global  nf;

% parameters for the SPG method %
gamma = 1e-4;            M = 2;
alphamin = 1e-8;         alphamax = 1;
f_1 = -inf*gpuArray(ones(M,1));    AtAb = A(A(b,1),2);

[f, AtAx] = func(A,b,mu,x);
g = grad(A,Atype,AtAb,mu,AtAx);
nf  = nf + 1;     iter_in = 0;

err = norm(soft_thresh(x-g,1) - x);
alphas = 1;

f_1(1) = f;
k = 1;
while (err/max(abs(f),1)) >= eps
  
%   if (k > maxiter)
%     break
%   end

  f_max = max(f_1);
  d = soft_thresh(x-alphas*g, alphas) - x;
  delta = sum(d.*g) + sum(abs(x+d)) - sum(abs(x));
  
  [x_new,f_new,g_new] = linesearch(A,Atype,b,AtAb,mu,f_max,gamma,delta,x,d);  
  f_1(mod(k,M)+1) = f_new;
  
  dx = x_new - x;
  dg = g_new - g;
  xdotg = sum(dx.*dg);
  xsqr =  norm(dx)^2;
  alphas = max(alphamin,min(alphamax,xsqr/xdotg));
  
  x = x_new;
  g = g_new;
  f = f_new;
  
  err = norm(soft_thresh(x-g,1) - x);
%   fprintf(' Iter_in = %3.0d    fval = %4.2f     err = %4.2f\n', iter_in, f, err)
  iter_in = iter_in + 1;
  k = k + 1;
  
end

% fprintf(' iter_in = %3.0d', k)


