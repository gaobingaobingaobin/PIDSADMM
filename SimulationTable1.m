% Total variation denoising with random data
%% Generate problem data

rand('seed', 0);
randn('seed', 0);

% n = 100;
       % number of examples
nn = 100;       % number of features
% pp = 100/n;      % sparsity density
dimension=10;
for ii=1:dimension
    n=nn+(ii-1)*100;
    nnn(ii)=n;
x0 = ones(n,1);
for j = 1:3
    idx = randsample(n,1);
    k = randsample(1:10,1);
    x0(ceil(idx/2):idx) = k*x0(ceil(idx/2):idx);
end
b = x0 + randn(n,1);

lambda = 5;

e = ones(n,1);
D = spdiags([e -e], 0:1, n,n);
DTD=D'*D;
[v,d]=eigs(DTD);
r=max(d(:));


  alpha=-0.1;
    % [x1 history1] = lassoLADMM(A, b, lambda, 1.0,r);
    [x3 history3] = total_variationILSADMM(b, lambda, 1.0,r,alpha);
    [x2 history2] = total_variationLSADMM(b, lambda, 1.0,r,alpha);
    % kk1(ii)=history1.iteration;
    % tt1(ii)=history1.time;
    kk2(ii)=history2.iteration;
    tt2(ii)=history2.time;
    kk3(ii)=history3.iteration;
    tt3(ii)=history3.time;

    alpha=0.1;
    % [x11 history11] = lassoLADMM(A, b, lambda, 1.0,r);
    [x33 history33] = total_variationILSADMM(b, lambda, 1.0,r,alpha);
    [x22 history22] = total_variationLSADMM(b, lambda, 1.0,r,alpha);
    % kk11(ii)=history11.iteration;
    % tt11(ii)=history11.time;
    kk22(ii)=history22.iteration;
    tt22(ii)=history22.time;
    kk33(ii)=history33.iteration;
    tt33(ii)=history33.time;
    fprintf('$  %3d$ & %3d(%10.2f)&%3d(%10.2f)&%10.2f(%10.2f)&%3d(%10.2f)&%3d(%10.2f)&%10.2f(%10.2f)\\\\ \n',...
        nnn(ii),...   % n
        kk2(ii), tt2(ii),kk3(ii), tt3(ii), kk3(ii)/kk2(ii),tt3(ii)/tt2(ii),...
        kk22(ii), tt22(ii),kk33(ii), tt33(ii), kk33(ii)/kk22(ii),tt33(ii)/tt22(ii));
    clear D
    clear DTD
    clear b
end