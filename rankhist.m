clear all
close all

rall=zeros(1000,3);
load toyresults0p2
figure(2),plot(sort([lambda_true(2:end) 0],'descend'),'b','LineWidth',2)
% hold on;
% plot(sort(lambda,'descend'),'g')
rall(:,1)=r(1001:end)';
% figure(1),hist(rshow)
% set(gca,'FontSize',16,'fontweight','b')
% set(gca, 'xtick', min(rshow):1:max(rshow))

load toyresults0p5
% hold on;
% plot(sort(lambda,'descend'),'b')
rall(:,2)=r(1001:end)';
% figure(2),hist(rshow)
% set(gca,'FontSize',16,'fontweight','b')
% set(gca, 'xtick', min(rshow):1:max(rshow))

load toyresults0p8
hold on;
plot(sort(lambda,'descend'),'r','LineWidth',2)
ylabel('\lambda','FontSize',20,'fontweight','b','rot', 0)
set(gca,'FontSize',16,'fontweight','b')
% legend('truth','20% training','50% training','80% training')
legend('truth','inferred')

rall(:,3)=r(1001:end)';
rall=rall+1.2;
% [y, b] = hist(rall);
% rshow=r(1001:end);
% figure(3),hist(rshow)
figure(1),hist(rall);
xlabel('Rank','FontSize',16,'fontweight','b')
set(gca,'FontSize',16,'fontweight','b')
set(gca, 'xtick', 17:1:23)
legend('20% training','50% training','80% training')
