function [U lambda pr llikevec time_trace U_all] = PTF_DistriParaOnlineGibbs(xi,id,R,iterN,isbatch,M,cutmode,trainFrac,miniFrac)
% cutmode: the mode along with we seperate the data
% M: number of machines used   
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

xitest=xi(idall(Train+1:end));
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
    subdimN{m}(cutmode)=idstartend(m,2)-idstartend(m,1)+1;
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
a=5e-1;
dir_a=cell(1,K);
dir_a{cutmode}=cell(1,M);% store dir_a, on the cutmode, it has M cells to store dir_a on each machine 
U=cell(1,K);
for k=1:K
    if k~=cutmode
        dir_a{k}=a*N(k)*ones(N(k),R)./max(N);
        U{k} = sampleDirMat(a*ones(1,N(k)),R);
        U{k} = U{k}';
    end
end
c=1;epsi=1/R;pr_a = c*epsi;pr_b = c*(1-epsi);pr=betarnd(pr_a,pr_b);
gr=0.1;lambda_a = gr*ones(1,R);lambda_b = (pr/(1-pr))*ones(1,R);lambda=gamrnd(gr,pr/(1-pr),1,R);      
if isbatch
    miniFrac=1;
end
    
    
%cutmode=1;

mypool=parpool(M);

dir_am=cell(1,M);
Um=cell(1,M);
parfor m=1:M
    dir_am{m}=a*ones(subdimN{m}(cutmode),R);
    for r=1:R
        Um{m}(:,r)=sampleDirMat(dir_am{m}(:,r)',1)'/M;
    end
end

llikevec=zeros(iterN,1);
t0 = 0;

subidsel=cell(1,M);% store the entry inddex in each machine, along the cutmode-th mode, the indices in all machines all start from 1        
subxisel=cell(1,M);% store the entry value in each machine  
xsumrsum=cell(1,K);
parfor m=1:M
    subidsel{m}=zeros(floor(miniFrac*subN(m)),K);
    subxisel{m}=zeros(floor(miniFrac*subN(m)),1);
end

xsumr=cell(1,M);
xr=cell(1,M);

dirtmp=cell(1,M);
for m=1:M
    dirtmp{m}=cell(1,K);   
    for k=1:K
        if k==cutmode
            dirtmp{m}{k}=a*ones(subdimN{m}(k),R);
        else
            dirtmp{m}{k}=a/M*ones(N(k),R);               
        end
    end 
end
xr_old=cell(1,M);
xsumr_old=cell(1,M);
    
    for iter=1:iterN
        tic
        if isbatch
            gam_t = 1;
        else
            gam_t = (iter+t0)^(-0.5);
        end
        zetair=cell(1,M);
        parfor m=1:M
            idid=randperm(subN(m),floor(miniFrac*subN(m)));
            subxisel{m}=xisub{m}(idid);
            for k=1:K
                subidsel{m}(:,k)=idsub{m,k}(idid);
                if k==cutmode
                    if m>1
                        subidsel{m}(:,k)=subidsel{m}(:,k)-idstartend(m,1)+1;
                    end
                end
            end
            zetair{m}=computezetair_parallel(U,Um{m},subidsel{m},lambda,cutmode);
            xiselectr=repmat(subxisel{m},1,R).*zetair{m};
            if iter==1
                [xsumr{m},xr{m}]=tensorsum(xiselectr,subidsel{m},subdimN{m});
                xr{m}=xr{m}/miniFrac;
                for k=1:K
                    xsumr{m}{k}=xsumr{m}{k}/miniFrac;
                end         
            else
                [xsumr_temp,xr_temp]=tensorsum(xiselectr,subidsel{m},subdimN{m});
                xr{m} = (1-gam_t)*xr_old{m} + gam_t*xr_temp/miniFrac;
                for k=1:K
                    xsumr{m}{k}=(1-gam_t)*xsumr_old{m}{k} + gam_t*xsumr_temp{k}/miniFrac;
                end              
            end
            dirtmp{m}=dirpara(xsumr{m},subidsel{m},a,cutmode,m,dirtmp{m});           
        end
        xrall=0;
        for m=1:M
            xrall=xrall+xr{m};
        end
        for k=1:K
            xsumrsum{k}=0;              
        end
        for k=1:K
            if k==cutmode
                for m=1:M
                    xsumrsum{k}=xsumrsum{k}+sum(dirtmp{m}{k});
                end                
            else
                for m=1:M
                    xsumrsum{k}=xsumrsum{k}+dirtmp{m}{k};
                end
            end
        end
        
        % update pr
        pr_a = c*epsi+xrall;
        pr_b = c*(1-epsi)+gr;
        pr = pr_a./(pr_a+pr_b);     

        % update lambda
        lambda_a = gr+xrall;
        lambda_b = pr;
        lambda = lambda_a.*lambda_b;
        
        for k=1:K
            if k==cutmode
                parfor m=1:M
                    dir_am{m}=dirtmp{m}{k};
                    Um{m}=dirtmp{m}{k}./repmat(xsumrsum{k},size(dir_am{m},1),1);
                end
            else
                dir_a{k}=xsumrsum{k};
                U{k}=dir_a{k}./repmat(sum(dir_a{k}),N(k),1);
            end
        end
        
        U{cutmode}=[];
        for m=1:M
            U{cutmode}=[U{cutmode};Um{m}];
        end
        
        xr_old = xr;
        xsumr_old = xsumr;

        if iter==1 
            time_trace(iter) = toc;
            tic;
        else
            time_trace(iter) = time_trace(iter-1) + toc;
            tic;
        end
        U_all=U;       
        [llike mae rmse mse]=evaluation(xitest,idtest,U_all,lambda);
        llikevec(iter)=llike;
        rmsevec(iter)=rmse;
        maevec(iter)=mae;
        msevec(iter)=mse;
        fprintf('iteration= %d;loglikelihood= %f, mae=%f, mse=%f, rmse=%f, time elapsed= %f\n', iter, llike, mae, mse, rmse, time_trace(iter));
        figure(1);
        subplot(2,1,1); plot(time_trace(1:iter),llikevec(1:iter));
        xlabel('Time (seconds)');
        ylabel('Heldout log-likelihood'); 
        subplot(2,1,2);plot(sort(lambda,'descend'));
        xlabel('Weights of rank-1 components');
        %figure(2),bar(xsumrsum{1})
        drawnow;
    end
    
    delete(mypool);
end