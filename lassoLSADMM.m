function [z, historz] = lassoLSADMM(A, b, lambda, rho,r,alpha)
% tau=(alpha^2-alpha+4)/(alpha^2-2*alpha+5);
% r=r*tau;
%r>rho||A^TA||
% lasso  Solve lasso problem via ADMM
%
% [z, historz] = lasso(A, b, lambda, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| Az - b ||_2^2 + \lambda ||z ||_1
%
% let x-Az=0;
% Langrange(x,z)=1/2*|| x - b ||_2^2 + \lambda ||z||_1 - u^T(x-Az) + rho/2*|| x - Az ||_2^2
% solution of x-subproblem
% x=1/(1+rho)*(b+u+rho*A*z)
%solution of z-subproblem
%z=S_{u/(tau*r)}(x+q/(tau*r) )where
%r>rho||A^TA|| and tau\in[0.8,1) and q=A'*（u-rho*(x-A*z)）
% The solution is returned in the vector z.
%
% historz is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (tzpical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~bozd/papers/distr_opt_stat_learning_admm.html
%

t_start = tic;
%% Global constants and defaults

QUIET    = 0;
MAX_ITER = 1000;
% MAX_ITER = 50;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
%% Data preprocessing

[m, n] = size(A);

% save a matrix-vector multiplz
Atb = A'*b;
%% ADMM solver

x = zeros(m,1);
z = zeros(n,1);
u = zeros(m,1);

% cache the factorization
% [L U] = factor(A, rho);

% if ~QUIET
%     fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
%         'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
% end
% ATA=A'*A;
% [v,d]=eigs(ATA);
% r=max(d(:));
% save r
for k = 1:MAX_ITER
    
    % x-update
    %     q = Atb + rho*(z - u);    % temporarz value
    %     if( m >= n )    % if skinnz
    %         x = U \ (L \ q);
    %     else            % if fat
    %         x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
    %     end
    % x=1/(1+rho)*(b+u-rho*A*z);错过无
    x=1/(1+rho)*(b+u+rho*A*z);
    % z-update with relaxation
    zold = z;
    %     x_hat = alpha*x + (1 - alpha)*zold;
    x_hat=x;
    u = u - alpha* rho*(x - A*z);% SADMM中第一个multiplier的parameter等于alpha
    %     z = shrinkage(x_hat + u, lambda/rho);
    %z=S_{u/(tau*r)}(x+q/(tau*r) )where
    %r>rho||A^TA|| and tau\in[0.8,1) and
    q=-A' * (u-rho*(x-A*z));
    z = shrinkage(z + q/r, lambda / r);
    % u-update
    u = u - rho*(x - A*z);
    
    % diagnostics, reporting, termination checks
    historz.objval(k)  = objective(A, b, lambda, x, z);
    
 historz.r_norm(k)  = norm(x - A*z);
    historz.s_norm(k)  = norm(-rho*A*(z - zold));
    
    historz.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-A*z));
    historz.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
    
%     if ~QUIET
%         fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
%             historz.r_norm(k), historz.eps_pri(k), ...
%             historz.s_norm(k), historz.eps_dual(k), historz.objval(k));
%     end
    
    if (historz.r_norm(k) < historz.eps_pri(k) && ...
            historz.s_norm(k) < historz.eps_dual(k))
             historz.iteration=k;
        historz.time=toc(t_start);
        break;
    end
    
end

% if ~QUIET
%     toc(t_start);
% end
end

function p = objective(A, b, lambda, x, z)
% p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) );
p = ( 1/2*sum((x - b).^2) + lambda*norm(z,1) );
end

function z = shrinkage(x, kappa)
z = max( 0, x - kappa ) - max( 0, -x - kappa );
end

function [L U] = factor(A, rho)
[m, n] = size(A);
if ( m >= n )    % if skinnz
    L = chol( A'*A + rho*speze(n), 'lower' );
else            % if fat
    L = chol( speze(m) + 1/rho*(A*A'), 'lower' );
end

% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');
end
