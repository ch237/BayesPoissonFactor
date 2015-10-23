function zetair = computezetair_parallel(U,Um,id,lambda,cutmode)%computezetair_parallel(U,id,lambda)
K=length(U);
% N=zeros(1,K);
R=size(Um,2);
% for k=1:K
%     if k==cutmode
%         
%     else
%         N(k)=size(U{1,k},1);
%     end
% end
if iscell(id)
else
    id_row=size(id,1);
    id=mat2cell(id,id_row,ones(1,K));
end
Nnon0=length(id{1,1});
log_zetair=repmat(log(lambda),Nnon0,1);
for k =1:K
    if k==cutmode
        log_zetair=log_zetair + log(Um(id{1,k},:));
    else
        log_zetair=log_zetair + log(U{1,k}(id{1,k},:));
    end
end
zetair=exp(log_zetair);

zetair=zetair./repmat(sum(zetair,2),1,R);
end