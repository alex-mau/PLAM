% Postprocess the solution

function xp = PostPro(x,X,y,eps)
%    x:  is the solution obtained by CDRSM;
%    X:  is the design matrix;
%    y:  is the observation vector;
%  eps:  is the tolerence of soloution, usually, it is 0.01;

p = size(X,2);
J = find(abs(x) >= eps);
xp = zeros(p,1);
AA = X(:,J);
bb = AA'*y;
xx = AA'*AA\bb;
xp(J) = xx;


%%  === Alternative processing =========
% b = X'*y;
% p = size(X,2);
% J = find(abs(x) >= eps);
% xp = zeros(p,1);
% AA = X(:,J);
% bb = b(J);
% xx = AA'*AA\bb;
% xp(J) = xx;

end