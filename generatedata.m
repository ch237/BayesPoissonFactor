clear all
xr=10*randperm(500,20);
xr=[xr,abs(randn(1,30))];
a=abs(2+randn(1,50));
lambda=gamrnd(xr,a);
K=3;
R=50;
N=[300,300,300];
Ut=cell(1,K);
for k=1:K
    Ut{k}=gamrnd(1e-1+1e-1*abs(randn(N(k),R)),1e-2);
    Ut{k}=Ut{k}./repmat(sum(Ut{k}),N(k),1);
end
xmean=zeros(N);
tensor=zeros(N);
for i=1:N(1)
    for j=1:N(2)
        for k=1:N(3)
            xmean(i,j,k)=sum(lambda.*Ut{1}(i,:).*Ut{2}(j,:).*Ut{3}(k,:));
%             tensor(i,j,k)=poissrnd(xmean(i,j,k));
        end
    end
end
disp('rate generated')
tensor=poissrnd(xmean);
clear xmean;
% for i=1:N(1)
%     for j=1:N(2)
%         for k=1:N(3)
%              tensor(i,j,k)=poissrnd(xmean(i,j,k));
%         end
%     end
% end
plot(lambda);
[r,c,v] = ind2sub(size(tensor),find(tensor> 0));
xi=zeros(length(r),1);
xi=zeros(length(r),1);
[rval rid]=sort(r,'ascend');
r=rval;
c=c(rid);
v=v(rid);
for i=1:length(r)
    xi(i)=tensor(r(i),c(i),v(i));
end
id=cell(1,K);
id{1}=r;
id{2}=c;
id{3}=v;
clear r
clear c
clear v

co_xr=10*randperm(500,10);
co_xr=[co_xr,abs(randn(1,40))];
co_a=abs(2+randn(1,50));
co_lambda=gamrnd(co_xr,co_a);
co_xmean=zeros(N(1),N(1));
for i=1:N(1)
    for j=1:N(1)
        co_xmean(i,j)=sum(co_lambda.*Ut{1}(i,:).*Ut{1}(j,:));
%             tensor(i,j,k)=poissrnd(xmean(i,j,k));
    end
end
bimat=poissrnd(co_xmean);
% bimat=((bimat>0)+(bimat'>0))>1;
bimat=((bimat>0)+(bimat'>0))>0;

[r,c] = ind2sub(size(bimat),find(bimat> 0));
co_xi=zeros(length(r),1);
[rval rid]=sort(r,'ascend');
r=rval;
c=c(rid);
for i=1:length(r)
    co_xi(i)=bimat(r(i),c(i));
end
co_id=cell(1,2);
co_id{1}=r;
co_id{2}=c;

clear r
clear c
save bimattensortoyall
save bimattensortoy co_id co_xi id xi