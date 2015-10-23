function dir_a=dirpara(xsumr,id,a,cutmode,M,dir_a)
K=length(xsumr);
% dir_a=xsumr;
% dir_a=cell(1,K);
% if iter==1   
%     for k=1:K
%         if k==cutmode
%         dir_a{k}=a*ones(size(xsumr{k}));
%         else
%             dir_a{k}=zeros(size(xsumr{k}));
%             if M==1
%                 dir_a{k}=a+dir_a{k};
%             end
%         end
%     end
% else
% %     dir_a=xsumr;
% %     for k=1:K
% %         dir_a{k}=a+
% %     end
% end
if iscell(id)
else
    id_row=size(id,1);
    id=mat2cell(id,id_row,ones(1,K));
end
R=size(xsumr{1},2);
for r=1:R
    for k=1:K
        if k==cutmode
            dir_a{k}(id{k},r) = a+xsumr{1,k}(id{k},r);
        else
%             if M==1
            dir_a{k}(id{k},r)=a/M+xsumr{1,k}(id{k},r);
%             end
        end
     end       
end

% if iter==1
% for r=1:R
%     for k=1:K
%         if k==cutmode
%             dir_a{k}(id{k},r) = dir_a{k}(id{k},r)+xsumr{1,k}(id{k},r);
%         else
% %             if M==1
%                 dir_a{k}(id{k},r)=dir_a{k}(id{k},r)+xsumr{1,k}(id{k},r);
% %             end
%         end
%      end       
% end
% else
%     for r=1:R
%         for k=1:K
%             if k==cutmode
%                 dir_a{k}(id{k},r) = a+xsumr{1,k}(id{k},r);
%             else
%     %             if M==1
%                     dir_a{k}(id{k},r)=a+xsumr{1,k}(id{k},r);
%     %             end
%             end
%          end       
%     end
% end

