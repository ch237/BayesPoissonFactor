clear all;close all;clc;
load xi;
load id;
len1=length(find(id{1,3}==1));
len2=length(find(id{1,3}==2394));
lenrmv=len1+len2;
len=length(xi)-lenrmv;
id_temp=cell(1,3);
for k=1:3
   id_temp{k}=zeros(lenrmv,1); 
end
xi_temp=zeros(lenrmv,1);
j=1;
for i=1:length(xi)
    if id{3}(i)~=1 && id{3}(i)~=2394
        xi_temp(j)=xi(i);
        for k=1:3
            id_temp{k}(j)=id{k}(i);
        end
        j=j+1;
    end
end
uaid=unique(id_temp{1});
[aid aname]=textread('/Users/changwei/Documents/Research/data/dukescholar/aidaname.txt','%s%s','delimiter',':');
aoidnid=zeros(length(aid),2);
j=1;
fid=fopen('aidanamenew.txt','w+');
for i=1:length(aid)
    aoidnid(i,1)=i;
    if any(uaid==i)
        aoidnid(i,2)=j;
        fprintf(fid,'%d:%s\r\n',j,aname{i});
        j=j+1;
    end
end
fclose(fid);

[vid vname]=textread('/Users/changwei/Documents/Research/data/dukescholar/vidvname.txt','%s%s','delimiter','|');
fid=fopen('vidvnamenew.txt','w+');
for i=1:length(vname)
    if i==1 | i==2394
    elseif i<2394
       fprintf(fid,'%d|%s\r\n',i-1,vname{i});
    elseif i>2394
       fprintf(fid,'%d|%s\r\n',i-2,vname{i}); 
    end
end
fclose(fid);
for i=1:length(xi_temp)
    temp=id_temp{1}(i);
    id_temp{1}(i)=aoidnid(temp,2);
    temp=id_temp{3}(i);
    if temp<2394
        id_temp{3}(i)=temp-1;
    else
        id_temp{3}(i)=temp-2;
    end
end

wordoldn=max(id{2});
wordocur=zeros(wordoldn,1);
worddococur=zeros(wordoldn,1);
for i=1:length(xi_temp)
    wordocur(id_temp{2}(i))=wordocur(id_temp{2}(i))+xi_temp(i);
    worddococur(id_temp{2}(i))=worddococur(id_temp{2}(i))+1;
end

WO=textread('vocab.txt','%s');
% clc
tf=zeros(size(wordocur));
tf(find(wordocur>0))=log(wordocur(find(wordocur>0)));
idf=zeros(size(wordocur));
idf(find(worddococur>0))=log(max(id_temp{1})*max(id_temp{3})./worddococur(find(worddococur>0)));
tfidf=tf.*idf;
% idf=wordocur./(worddococur+eps)
% [a b]=sort(tfidf,'descend');
% [a b]=sort(tfidf);
% [a b]=sort(idf,'descend');

% [a b]=sort(wordocur,'descend');
% wordocurtemp=wordocur;
% wordocurtemp(find(wordocur==0))=inf;
[a b]=sort(worddococur,'descend');
% [a b]=sort(wordocur);
% idf(find(worddococur==0))=inf;
% [a b]=sort(idf);
% wordrmv=cell(1,1);
rmvn=100;
% wordrmv=zeros(rmvn,1);
str=[];
for i=1:rmvn
%     wordrmv{1}(i)=WO{b(i)};
    if mod(i,15)==0 | i==rmvn
        disp(str)
        str=[];
    else
        str=[str, ' ', WO{b(i)}];
    end
end


woidnid=zeros(max(id{2}),2);



mannualid=[7251 1298 7253 3003 3004 3706 5264 7450];
for i=1:length(mannualid)
    worddococur(mannualid(i))=0;
end

worddococur_temp=worddococur;
worddococur_temp(b(1:rmvn))=0;
j=1;
fid=fopen('vocabnew.txt','w+');
for i=1:length(worddococur_temp)
    woidnid(i,1)=i;
    if worddococur_temp(i)>0
        woidnid(i,2)=j;
        fprintf(fid,'%s\r\n',WO{i});
        j=j+1;
    end
end
fclose(fid);
for i=1:length(xi_temp)
    if woidnid(id_temp{2}(i),2)==0
        xi_temp(i)=0;
    end
end
zeron=length(find(xi_temp==0));
xi=zeros(length(xi_temp)-zeron,1);
id=cell(1,3);
for k=1:3
    id{k}=zeros(size(xi));
end
j=1;
for i=1:length(xi_temp)
    if xi_temp(i)==0
    else
        xi(j)=xi_temp(i);
        id{1}(j)=id_temp{1}(i);
        id{3}(j)=id_temp{3}(i);
        id{2}(j)=woidnid(id_temp{2}(i),2);
        j=j+1;
    end        
end
save xiidnew xi id

% subplot(3,1,1),bar(wordocur)
% subplot(3,1,2),bar(worddococur)
% subplot(3,1,3),bar(idf)