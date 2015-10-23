% function zetair = VBcomputezetair(xsumr,id,gr,xr,pr,a)
% K=length(xsumr);
% N=zeros(1,K);
% R=size(xsumr{1,1},2);
% for k=1:K
%     N(k)=size(xsumr{1,k},1);
% end
% Nnon0=length(id{1,1});
% zetair=ones(Nnon0,R);
% 
% for k =1:K    
%     expdirpart2=exp(repmat(psi(sum(xsumr{1,k})+N(k)*a(k)),Nnon0,1));
%     zetair=zetair.*exp(psi(xsumr{1,k}(id{1,k},:)+a(k)))./expdirpart2;
% end
% 
% zetair=zetair.*repmat(exp(psi(gr+xr)).*pr,Nnon0,1);
% % zetair=zetair;
% zetair=zetair./repmat(sum(zetair,2),1,R);
% end



function zetair = VBcomputezetair1(xsumr,id,gr,xr,pr,a)
K=length(xsumr);
N=zeros(1,K);
R=size(xsumr{1,1},2);
for k=1:K
    N(k)=size(xsumr{1,k},1);
end
if iscell(id)
else
    idtemp=id;
    id=cell(1,K);
    for k=1:K
        id{k}=idtemp(:,k);
    end
end
Nnon0=length(id{1,1});
% zetair=ones(Nnon0,R);
zetair=zeros(Nnon0,R);

for k =1:K    
    zetair=zetair+psi(xsumr{1,k}(id{1,k},:)+a(k));
    xsumcolr=psi(sum(xsumr{1,k})+N(k)*a(k));
    zetair = zetair - repmat(xsumcolr,Nnon0,1);
    %expdirpart2=expdirpart2+psi(xsumcolr);
end
% for k =1:K    
%     zetair=zetair.*exp(psi(xsumr{1,k}(id{1,k},:)+a(k)));
% end
% expdirpart2=0; 
% for k=1:K
%     expdirpart2=expdirpart2+psi(sum(xsumr{1,k})+N(k)*a(k));
% end
zetair = zetair + repmat(psi(gr+xr)+log(pr),Nnon0,1);
%zetair=exp(zetair+repmat(expgampart-expdirpart2,Nnon0,1));
% zetair=zetair+repmat(expgampart-expdirpart2,Nnon0,1);
% zetair=zetair.*repmat(exp(psi(gr+xr)).*pr,Nnon0,1);
% zetair=zetair;
zetair = exp(zetair);
zetair=zetair./repmat(sum(zetair,2),1,R);
end