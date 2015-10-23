function xi=xigenerate(X,id)
    K=length(id);
    Nnon0=length(id{1,1});
    xi=zeros(Nnon0,1);
    tempid=zeros(1,K);
    for i=1:Nnon0
        for k=1:K
            tempid(k)=id{1,k}(i);
        end
        if K==2
            xi(i)=X(tempid(1),tempid(2));
        elseif K==3
            xi(i)=X(tempid(1),tempid(2),tempid(3));
        end
    end
end