function [xsumr,xr]=tensorsum(xir,id,N)
[Nnon0,R]=size(xir);
K=length(N);
if iscell(id)
else
    id_row=size(id,1);
    id=mat2cell(id,id_row,ones(1,K));
end
xsumr=cell(1,K);

for k=1:K
    xsumr{1,k}=zeros(N(k),R);
end

idtemp=zeros(1,K);

for i=1:Nnon0
    for k=1:K
        idtemp(k)=id{1,k}(i);
    end
    for k=1:K
        xsumr{1,k}(idtemp(k),:)=xsumr{1,k}(idtemp(k),:)+xir(i,:);
    end
end

xr=sum(xir,1);

