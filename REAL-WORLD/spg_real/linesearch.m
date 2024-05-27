% nonmonotone line search used in the SPG method
function [x_new,f_new,g_new] = linesearch(A,Atype,b,AtAb,mu,f_max,gamma,delta,x,d)

global nf;

alphas = 1;

while (1)
  
  x_new = x + alphas*d;
  [f_new, AtAx] = func(A,b,mu,x_new);
  nf  = nf + 1;
  if (f_new <= f_max + gamma*alphas*delta) || (alphas <= 1e-8)
    break
  else
    alphas = alphas/2;
  end
  
end
g_new = grad(A,Atype,AtAb,mu,AtAx);
