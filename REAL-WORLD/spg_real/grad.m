function g = grad(A,Atype,AtAb,mu,AtAx)

%Calculate Gradient
if (Atype == 'I')
  g = mu*(AtAx - AtAb);
else
  g = mu*(A(A(AtAx,1),2) - AtAb);
end