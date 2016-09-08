% L1-regularized least-squares example
%% Generate problem data

randn('seed', 0);
rand('seed',0);

mm = 150;       % number of examples
nn = 500;       % number of features
% pp = 100/n;      % sparsity density
dimension=5;
for ii=1:dimension
    m=mm+(ii-1)*150+750;
    mmm(ii)=m;
    n=nn+(ii-1)*500+2500;
    nnn(ii)=n;
    p=100/n;
    x0 = sprandn(n,1,p);
    A = randn(m,n);
    A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
    % % save A
    % load A
    b = A*x0 + sqrt(0.001)*randn(m,1);
    % save b
    % load b
    ATA=A'*A;
    [v,d]=eigs(ATA);
    r=max(d(:));
    % save r %大约是7.9
    % load r
    lambda_max = norm( A'*b, 'inf' );
    lambda = 0.1*lambda_max;
    %% Solve problem
    % tau=0.8;
    alpha=-0.3;
    % [x1 history1] = lassoLADMM(A, b, lambda, 1.0,r);
    [x3 history3] = lassoILSADMM(A, b, lambda, 1.0,r,alpha);
    [x2 history2] = lassoLSADMM(A, b, lambda, 1.0,r,alpha);
    % kk1(ii)=history1.iteration;
    % tt1(ii)=history1.time;
    kk2(ii)=history2.iteration;
    tt2(ii)=history2.time;
    kk3(ii)=history3.iteration;
    tt3(ii)=history3.time;

    alpha=0.3;
    % [x11 history11] = lassoLADMM(A, b, lambda, 1.0,r);
    [x33 history33] = lassoILSADMM(A, b, lambda, 1.0,r,alpha);
    [x22 history22] = lassoLSADMM(A, b, lambda, 1.0,r,alpha);
    % kk11(ii)=history11.iteration;
    % tt11(ii)=history11.time;
    kk22(ii)=history22.iteration;
    tt22(ii)=history22.time;
    kk33(ii)=history33.iteration;
    tt33(ii)=history33.time;
    fprintf('$ %3d \\times %3d$ & %3d(%10.2f)&%3d(%10.2f)&%10.2f(%10.2f)&%3d(%10.2f)&%3d(%10.2f)&%10.2f(%10.2f)\\\\ \n',...
        mmm(ii),nnn(ii),...   %m * n
        kk2(ii), tt2(ii),kk3(ii), tt3(ii), kk3(ii)/kk2(ii),tt3(ii)/tt2(ii),...
        kk22(ii), tt22(ii),kk33(ii), tt33(ii), kk33(ii)/kk22(ii),tt33(ii)/tt22(ii));
    clear A 
    clear b
end
