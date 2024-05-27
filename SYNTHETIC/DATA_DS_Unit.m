% Generate the data for Dantzig selector
% Refer to: Dissussion: The Dantzig selector: statistical estimation when p
% is much larger than n;
% The Annals of Statistics, 2007, 6: 2385--2391;
% M. Friedlander and M. A. Saunders

% min  || \beta ||_1
% s.t. || X^T (y-X\beta) ||_{\infty} \leq \delta

function [X,D,y,betas,delta] = DATA_DS_Unit(n,p,T,sigma,Type,Seed)

%% === Wang code ===== Unit coloumn ====
% random generate X and normalize the column

if strcmp(Seed,'fix')
    rand('state',0);  randn('state',0);
end

X = randn(n,p);
for k = 1:p
    X(:,k) = X(:,k)/norm(X(:,k));
end

q = randperm(p);
q = q(1:T);
betas = zeros(p,1);
if strcmp(Type, 'I')
    betas(q) = randsrc(T,1).*( 1 + abs(randn(T,1)) ) ;    %% Case I
else
    betas(q) = (2*rand(T,1)-1).*( 1+abs(randn(T,1)) ) ; %% Case II
end

noise = sigma.*randn(n,1);
y = X*betas + noise;
delta = sqrt(2*log(p))*sigma;
D = ones(p,1);

clear noise perm  k q

% %%%======================================================
% randn('seed', 0);    rand('state',0);
%
% % sigma = 0.01;
% q = randperm(p);
% q = q(1:T);
% betas = zeros(p,1);
% betas(q) = (2*rand(T,1)-1).*( 1+abs(randn(T,1)) ) ;
% betas = sparse(betas);
%
% % betas = sprand(p,1,0.01);
% B = randn(n,p);
% X = B*spdiags(1./sqrt(sum(B.^2))',0,p,p); % normalize columns
% noise = sigma*randn(n,1);
% y = X*betas + noise;
% delta = sigma*sqrt(2*log(p));
% % D = sqrt( sum( X.^2, 1) )';
% D = ones(p,1);
% clear q B noise