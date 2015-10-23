function [llike mae rmse mse]=evaluation(xi,id_temp,U,lambda)
K=length(U);
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
% llike=sum(log(poisspdf(xi,zetai)+eps));
% llike=sum(-zetai+xi.*log(zetai)-log(factorial(xi)));
llike=sum(-zetai+xi.*log(zetai)-gammaln(xi+1));
rmse=norm(xi-zetai,'fro')/(sqrt(size(xi,1) * size(xi,2)));
mae=sum(abs(xi-zetai))/Np;
mse=rmse^2;

% errornorm2=sum((xi-zetai).^2)/Np;%sum(xi.^2);
