function r = sampleDirMat(a,n)
%a:input row vector, n: number of dirichlet samples
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);