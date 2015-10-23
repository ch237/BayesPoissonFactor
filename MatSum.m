function [Xpk,Xnk,Xk]=MatSum(P,N,K,id1,id2,Xpik)
Xpk=zeros(P,K);
Xnk=zeros(N,K);
% Xk=zeros(1,K);
for i = 1:length(id1)
    Xpk(id1(i),:)=Xpk(id1(i),:)+Xpik(i,:);
    Xnk(id2(i),:)=Xnk(id2(i),:)+Xpik(i,:);
end
Xk=sum(Xpk);
Xk=Xk';
