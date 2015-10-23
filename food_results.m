clear;
load('food_res_1');
store_info = csvread('3D_storedict.csv');
store_info(:,2) = store_info(:,2) + 1;
[a b] = sort(store_info(:,2));
store_info(:,1) = store_info(b,1);
store_info(:,2) = a;
item_info = csvread('3D_itemdict.csv');
item_info(:,2) = item_info(:,2) + 1;
[a b] = sort(item_info(:,2));
item_info(:,1) = item_info(b,1);
item_info(:,2) = a;

R = size(U{2},2);
num_top = 200;
for r=1:R
    [a b] = sort(U{2}(:,r),'descend');
    store_clusters(r,:) = store_info(b(1:num_top),1);
    [a b] = sort(U{3}(:,r),'descend');
    item_clusters(r,:) = item_info(b(1:num_top),1);    
end

