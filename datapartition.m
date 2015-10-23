function datapartition(xi,id,M,cutmode,trainFrac)

K=length(id);
for k=1:K 
    N(k) = max(id{k}); 
end 
Nnon0=length(xi);
Train=round(trainFrac*Nnon0);
idall=randperm(Nnon0);
idtest=cell(1,K);
idtrain=cell(1,K);
for k=1:K
    idtest{k} = id{k}(idall(Train+1:end));
    idtrain{k}=id{k}(idall(1:Train));
end

% xitest=xi(idall(Train+1:end));
xitrain=xi(idall(1:Train));  

xisub=cell(1,M);% store tensor entry values in each machine
idsub=cell(M,K);% store id of entries in each machine
idstartend=zeros(M,2);%store the end id for each machine
idstartend(1,1)=1;
idstartend(1,2)=floor(N(cutmode)/M);
for m=2:M
    idstartend(m,1)=idstartend(m-1,2)+1;
    if m<M
        idstartend(m,2)=floor(N(cutmode)*m/M);
    else
        idstartend(m,2)=N(cutmode);
    end
end
subdimN=cell(1,M);% store the dimenstionality of tensors in each machine
for m=1:M
    subdimN{m}=N;
%     if m==1
%         subdimN{m}(cutmode)=idstartend(m,2);
%     else
        subdimN{m}(cutmode)=idstartend(m,2)-idstartend(m,1)+1;
%     end
end
subN=zeros(M,1);% number of nonzero entries in each machine
for i=1:length(xitrain)
    if mod(i,1e5)==0
        fprintf('%d-th nonzero entry\n',i);
    end
    bin_id=find(idtrain{cutmode}(i)<=idstartend(:,2),1);
    subN(bin_id)=subN(bin_id)+1;
end
for m=1:M
    for k=1:K
        idsub{m,k}=zeros(subN(m),1);
    end
    xisub{m}=zeros(subN(m),1);      
end

subN=zeros(M,1);
for i=1:length(xitrain)
    bin_id=find(idtrain{cutmode}(i)<=idstartend(:,2),1);
    if mod(i,1e5)==0
        fprintf('Assign %d-th nonzero entry to machine %d\n',i,bin_id);
    end        
    subN(bin_id)=subN(bin_id)+1;
    for k=1:K
        idsub{bin_id,k}(subN(bin_id))=idtrain{k}(i);
    end
    xisub{bin_id}(subN(bin_id))=xitrain(i);
end
tempxisub=xisub;
tempidsub=idsub;
tempsubdimN=subdimN;
tempidstartend=idstartend;
idsub=cell(1,K);
for m=1:M
    xisub=tempxisub{m};
    for k=1:K
        idsub{k}=tempidsub{m,k};
    end
    subdimN=tempsubdimN{m};
    idstartend=tempidstartend(m,:);
    eval(sprintf('save partition%d xisub idsub subdimN idstartend',m))
end