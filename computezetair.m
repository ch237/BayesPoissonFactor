function zetair = computezetair(U,id,lambda)
K=length(U);
N=zeros(1,K);
R=size(U{1,1},2);
for k=1:K
    N(k)=size(U{1,k},1);
end
if iscell(id)
else
    id_row=size(id,1);
    id=mat2cell(id,id_row,ones(1,K));
end
Nnon0=length(id{1,1});
zetair=ones(Nnon0,R).*repmat(lambda,Nnon0,1);
for k =1:K
    zetair=zetair.*U{1,k}(id{1,k},:);
end
zetair=zetair+eps;
zetair=zetair./repmat(sum(zetair,2),1,R);
end