%Algorithm for solution to tiff ODEs using orthgonal collaction on finite
%elements

%N,i - number of elements
%K,j or k - number of collocation points within element
%pol - polynomial name or type whose roots form the points
%h - size of each element

%inputs
%x0 - intial point of integration
%xf - final point of integration
%x - total integration horizon

%z

N = 1;
K = 3;
pol = 'Gauss-Radau';

% x = xf-x0;
% h = x/N;

for iel = 1:N
    for jcp = 1:K        
    end
end

z = fsolve(@func1,[0 0 0]);
%rhs matrix
% func = @(x)(x^2-2*x+1);
% x = [0;0.155051;0.644949;1];
% K = 3;
% L = zeros(K+1,1);
% Ldot = zeros(K+1,1);
% a = zeros(K+1,1);
% for jk = 1:K
%     for j = 0:K        
%         [L(j+1),Ldot(j+1)] = Lagrange(x(jk+1),j,x);
%         a(j+1) = Ldot(j+1);
%     end
%     A(jk,:) = a';
%     b(jk) = func(
% end




% L = (tau-x(2))/(x(1)-x(2)).(tau-x(3))/(x(1)-x(3))