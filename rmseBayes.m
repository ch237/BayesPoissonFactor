function e=rmseBayes(xi,id_temp,U,lambda)
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
e=norm(xi-zetai,'fro')/(sqrt(size(xi,1) * size(xi,2) -1));

% function e = rmse(x,y)
% %RMSE    Standard deviation.
% %   For vectors, RMSE(X,Y) returns the root mean square error of X with
% %   respect to the true value Y.
% %
% %   Matrices are treated as vectors.
% %   
% %   RMSE(X,Y) normalizes by (N-1) where N is the sequence length.  This
% %   makes RMSE(X,Y)^2 the best unbiased estimate of the variance if X
% %   is a sample from a normal distribution centered in Y.
% %
% 
% % note:  typical residual in many problems is  norm(x-y,'fro')^2
% 
% % Algorithm ref.: none
% %
% % Author: A. Fusiello, 2006
% 
% 
% % The Frobenius norm is the square root of the sum of squares of the
% % entries. The typical residual in many problems is norm(x1-x2,'fro')^2
% 
% e = norm(x-y,'fro')/(sqrt(size(x,1) * size(x,2) -1));
% 
% % for vectors this is equivalent to:
% % e = sqrt(sum((x - y).^2)/(length(x)-1));