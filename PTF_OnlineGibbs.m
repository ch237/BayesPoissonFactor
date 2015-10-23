function [U lambda pr llikevec time_trace] = PTF_OnlineGibbs(xi,id,R,batchsize,numiters,isbatch)
    K=length(id);
    for k=1:K 
        N(k) = max(id{k}); 
    end
    
    % initial parameters of the dirichlet on U
    a=0.5*(min(N))*(ones(1,K)./N);
    for k=1:K
        for r=1:R
            dir_a{k,r} = a(k)*ones(N(k),1);
        end
    end    
    U=cell(1,K);
    for k=1:K
        U{1,k} = sampleDirMat(a(k)*ones(1,N(k)),R);
        U{1,k} = U{1,k}';
    end
    
    % initial parameters of the beta on pr
    c=1;
    epsi=1/R;    
    pr_a = c*epsi;
    pr_b = c*(1-epsi);
    pr=betarnd(pr_a,pr_b);
    
    % initial parameters of the gamma on lambda
    gr=0.1;
    lambda_a = gr*ones(1,R);
    lambda_b = (pr/(1-pr))*ones(1,R);
    lambda=gamrnd(gr,pr/(1-pr),1,R); 
    
    Nnon0=length(id{1,1});
    Train=round(0.95*Nnon0);
    idall=randperm(Nnon0);
    idtest=cell(1,K);
    for k=1:K
        idtest{k} = id{k}(idall(Train+1:end)); 
    end  
    
    perc=1;% fraction of total training data
    
    iterN=numiters;
    if isbatch
        Np = Train;
    else
        Np=batchsize;
%         Np=floor(Train*perc*0.1);
    end
    llikevec=zeros(iterN,1);
    t0 = 0;
    tic;
    for iter=1:iterN
        if isbatch
            gam_t = 1;
        else
            gam_t = (iter+t0)^(-0.5);
        end
        
        %if iter==1
        idselect = [];
%         n_tmp = randperm(Train);
        n_tmp = randperm(floor(Train*perc));
        idid = n_tmp(1:Np);
        idid=idall(idid);

        xiselect = xi(idid);
        for k=1:K
            idselect(:, k) = id{1, k}(idid);      
        end
        %end

        % compute the mean of xir's
        if iter==1
            zetair=computezetair_new(U,idselect,lambda);
            xiselectr=repmat(xiselect,1,R).*zetair;
            [xsumr,xr]=tensorsum(xiselectr,idselect,N);
            % update the sufficient statistics using the xir's
            % note: multiplying by Train/Np is a hack but seems to give
            % better result and faster convergence
            xr=xr*Train/Np;
            for k=1:K
                xsumr{k}= xsumr{k}*Train/Np;
            end
        else
            zetair=computezetair_new(U,idselect,lambda);
            xiselectr=repmat(xiselect,1,R).*zetair;
            [xsumr_temp,xr_temp]=tensorsum(xiselectr,idselect,N);

            % update the sufficient statistics using the xir's
            % note: multiplying by Train/Np is a hack but seems to give
            % better result and faster convergence            
            xr = (1-gam_t)*xr_old + gam_t*xr_temp*Train/Np;
            for k=1:K
                xsumr{k}=(1-gam_t)*xsumr_old{k} + gam_t*xsumr_temp{k}*Train/Np;
            end
        end

        % update pr
        pr_a = c*epsi+xr;
        pr_b = c*(1-epsi)+gr;
        pr = pr_a./(pr_a+pr_b);     

        % update lambda
        lambda_a = gr+xr;
        lambda_b = pr;
        lambda = lambda_a.*lambda_b;
        
        % update U
        for r=1:R
            for k=1:K
                dir_a{k,r}(idselect(:, k)) = a(k)+xsumr{1,k}(idselect(:, k),r);
                U{1,k}(:,r) = dir_a{k,r}'/sum(dir_a{k,r}); 
                % TBD: U changed. Should we re-estimate the sufficient statistics?               
             end       
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
 
        [llike mae rmse mse]=evaluation(xi(idall(Train+1:end)),idtest,U,lambda);
        llikevec(iter)=llike;
        rmsevec(iter)=rmse;
        maevec(iter)=mae;
        msevec(iter)=mse;
        fprintf('iteration= %d;loglikelihood= %f, mae=%f, mse=%f, rmse=%f, time elapsed= %f\n', iter, llike, mae, mse, rmse, time_trace(iter));
           
        %llike=loglike(xi(idall(Train+1:end)),idtest,U,lambda);
        %llikevec(iter)=llike;
        subplot(2,1,1); plot(time_trace(1:iter),llikevec(1:iter));
        xlabel('Time (seconds)');
        ylabel('Heldout log-likelihood'); 
        subplot(2,1,2);plot(sort(lambda,'descend'));
        xlabel('Weights of rank-1 components');
        %ylabel('Heldout log-likelihood');        
        drawnow;
        %fprintf('iteration= %d;loglikelihood= %f, time elapsed= %f\n', iter, llike, time_trace(iter)); 
    end
end