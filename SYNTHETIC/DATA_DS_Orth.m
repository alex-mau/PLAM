% Generate the data for Dantzig selector
% Refer to: Dissussion: The Dantzig selector: statistical estimation when p
% is much larger than n;
% The Annals of Statistics, 2007, 6: 2385--2391;
% M. Friedlander and M. A. Saunders

% min  || \beta ||_1
% s.t. || X^T (y-X\beta) ||_{\infty} \leq \delta

function [X,D,y,betas,delta] = DATA_DS_Orth(n,p,T,sigma,Type,Seed)

% %%============== Wang Code == orthogonal rows=======
if strcmp(Seed,'fix')
    rand('state',0);  randn('state',0);
end

X = randn(n,p);   [Q,~] = qr(X');    X = Q(:,1:n)';

q = randperm(p);    q = q(1:T);
betas = zeros(p,1);
if strcmp(Type, 'I')
    betas(q) = randsrc(T,1).*( 1 + abs(randn(T,1)) ) ;    %% Case I
else
    betas(q) = (2*rand(T,1)-1).*( 1+abs(randn(T,1)) ) ; %% Case II
end

noise = sigma.*randn(n,1);
y = X*betas + noise;
delta = sqrt(2*log(p))*sigma;
D = zeros(p,1);
for ki = 1 : p
    D(ki) = norm( X(:,ki) );
end
clear R Q q noise ki

% % %%%%%=====================================================================
% % %
% rand('state',0);   randn('state',0);
%
% q = randperm(p);             % p is the dimension of \beta;
% q = q(1:T);                  % T is the number of nonzero components;
% betas = zeros(p,1);          % The given solution \beta;
% % betas(q) = sign(randn(T,1));
% betas(q) = randsrc(T,1).*( 1 + abs(randn(T,1)) ) ;
% [X,R] = qr(randn(p,n),0);
% X = X';                      % X*X' = I
% y = X*betas + sigma*randn(n,1);
% delta = sigma*sqrt(2*log(p));
% D = zeros(p,1);
% for ki = 1 : p
%     D(ki) = norm( X(:,ki) );
% end
% clear q R
