function [f, AtAx] = func(A,b,mu,x)

%Evaluate the function value
AtAx = A(A(x,1),2);
f = sum(abs(x)) + norm(AtAx - b)^2*mu/2;







