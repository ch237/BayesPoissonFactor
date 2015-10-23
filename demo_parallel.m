clear all;
load xi_id

R=50; % rank
iterN=200; % number of MCMC iterations
M=4;% number of splits of dataset, set it equal or less than the # of cores, set this parameter only when use the parallel implementation 
cutmode=2; % along which mode to split the data, set this parameter only when use the parallel implementation
isbatch=0; % whether online or batch Gibbs (0 is online) 
trainFrac=0.95; % fraction of data used to train the model
miniFrac=0.05; % fraction of training data used in each iteration of the online algorithm
paralellFlag=0;%1: use the paralel and online implementation; 0: use the online implementation

if paralellFlag==1
    [U lambda pr llikevec time_trace U_all] = PTF_DistriParaOnlineGibbs(xi,id,R,iterN,isbatch,M,cutmode,trainFrac,miniFrac);
    save results U lambda pr llikevec time_trace U_all
elseif paralellFlag==0
    batchsize=floor(length(xi)*trainFrac*miniFrac);
    [U lambda pr llikevec time_trace] = PTF_OnlineGibbs(xi,id,R,batchsize,iterN,isbatch);
    save results U lambda pr llikevec time_trace
end

xid=1;
while(xid~=0)
xid=input('Please input the topic you want show, and input 0 to exit:');
fprintf(' xid=%f\n ',xid);
end