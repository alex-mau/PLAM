% This codes solve the Dantzig selector problem, which can be equivalently formulated as:
%
%  min_x  ||x||_1  subject to  A'(Ax - b) = y and ||diag.*y||_\infty <= delta,
% 
% by alternating direction method (see Lu and Zhang (2009)). We solve the subproblem:
%
%  min_x ||x||_1 + (mu/2)*||A'(Ax - b)||_^2
%
% by the Nonmonotone Spectral Projected Gradient method for non-smooth function
% (see "An Augmented Lagrangian Approach for Sparse Principal Component Analysis"
% by Lu and Zhang).
%
% Call: [xsol, iter, gap, time] = selector(A,typeA,b,delta,eps,maxiter)
%
% input:
%
% A       - n x p matrix
% diag    - p x 1 positive vector for scaling.
% b       - n x 1 vector of observations.
% delta   - scalar input
% Atype   - type of the matrix A: ' ' or 'I'; 'I' means A*A' = I.
% eps     - tolerance for alternating direction method
% maxiter - maximum number of iterations
%
% output:
%
% xsol - approx solution to Dantzig selector problem
% gap  - primal-dual obj gap
% iter - number of iterations for the alogrithm
% time - CPU time for the algorithm
%
% Written by Zhaosong Lu, Ting Kei Pong and Yong Zhang
% Nov 13, 2010


function [x, iter, dval, time,Iter_inner] = selector_real(A,diag,Atype,b,delta,eps,maxiter)

% fprintf(' ************ start the iterations ************\n')

global iter_in
iter_in = 1;

tic;

[~,p] = size(A);

% Create function handle
if ~isa(A,'function_handle')
  A = @(x,mode) explicitMatrix(A,x,mode);
end

Atb = A(b,2);

% Initialization
iter = 1;  Iter_inner = 0;
x = zeros(p,1);
lambda = zeros(p,1);
mu = max(10/delta/sqrt(p), 200/sqrt(2*p*log(p)));
if Atype == 'I'
  mu = max(1/delta, 20/sqrt(2*log(p)));
end
AtAx = A(A(x,1), 2);     dval = -inf;

% Main loop
while iter <= maxiter
  
  % Solve the alternating direction for y
  tmp = lambda/mu;
  y = min(max(AtAx + tmp, -delta./diag+Atb),delta./diag+Atb);  
  
  % Solve the alternating direction for x using SPG
  eps2 = max(eps/100,100/iter^2);
  [x,Inner] = spg2(A,Atype,y-tmp,mu,x,eps2,maxiter);
  Iter_inner = Iter_inner + Inner;
  AtAx = A(A(x,1), 2);
  res = AtAx - y;
    
  % Calculating function values
  f_new = norm(x,1);
  
  tmp1 = norm(A(A(lambda,1),2),'inf');
  if tmp1 <= 1
    lambdabar = lambda;
  else lambdabar = lambda/tmp1;
  end
  dval = max(-delta*norm(lambdabar./diag,1)-lambdabar'*Atb, dval);
   
  pfeas = (norm(diag.*(AtAx - Atb),'inf')-delta)/max(norm(x,2),1);
    
  % Stopping criterion
  gap = abs(f_new - dval)/max(f_new,1);
  if max(gap/4,pfeas) <= eps
%     dfeas = (norm(A(A(lambdabar,1),2),'inf')-1)/max(norm(lambdabar,2),1);
%     fprintf(' Iter = %3.0d    fval = %7.3f   dval = %7.3f   pfeas = %6.4f   dfeas = %6.4f   time = %5.1f\n', iter, f_new, dval, pfeas, dfeas, toc)
    break;
  end
  
  % Update Lagrangian mutiplier
  lambda = lambda + mu*res;
%   dfeas = (norm(A(A(lambdabar,1),2),'inf')-1)/max(norm(lambdabar,2),1);
%   fprintf(' Iter = %3.0d    fval = %7.3f   dval = %7.3f   pfeas = %6.4f   dfeas = %6.4f   time = %5.1f\n', iter, f_new, dval, pfeas, dfeas, toc)
  iter = iter + 1;

end
time = toc;