function errornorm2=L2error(xi,id_temp,V,hr)
K=length(V);
U=cell(1,K);
for k=1:K
    U{1,k}=exp(V{1,k});
    Nk=size(U{1,k},1);
    U{1,k}=U{1,k}./repmat(sum(U{1,k}),Nk,1);
end
lambda=exp(hr);
if iscell(id_temp)
    id=id_temp;
else
id=cell(1,K);
for k=1:K
    id{1,k}=id_temp(:,k);
end
end


R=size(U{1,1},2);
Np=length(id{1,1});
zetair=ones(Np,R).*repmat(lambda,Np,1);
for k =1:K
    zetair=zetair.*U{1,k}(id{1,k},:);
end
zetai=sum(zetair,2);
errornorm2=sum((xi-zetai).^2)/Np;%sum(xi.^2);
