function [Bel,temp_Bel,Neighbor_num,Neighborhood]=DS(Neighbor,dataD,gamma_ds,W,distance)
%% 使用最大最小标准化对相似性矩阵W进行归一化
[~,b]=size(W);
for i=1:b
    for j=1:b
        if i==j
            Sim(i,j)=1;
        elseif i<j
            Sim(i,j)=W(i,j);
        else
            Sim(i,j)=W(i-1,j);
        end
    end
end
%Sim是每列形成一个样本的相似度，故归一化是按照列进行归一化
Sim=(mapminmax(Sim',0,1))';

%% DS 证据理论
w1=0.5;
w2=0.5;
value=0.95;

[N,~]=size(Neighbor);
u=unique(dataD);
[y,~]=size(u);
temp_Bel=zeros(N,y);
for i=1:N  %对于N个样本,每列作为一个样本的邻域集
    Neighborhood=find(Neighbor(:,i)==1);
    neighbor_ins=setdiff(Neighborhood,i);
    [d,~]=size(neighbor_ins);
    Neighbor_num(i)=d;
    
    if ~isempty(neighbor_ins)
        neighbor_label=dataD(neighbor_ins')';
        target_label=dataD(i);
        target_label_ID=find(gamma_ds(1,:)==target_label);
        target_gamma_ds = gamma_ds(2,gamma_ds(1,:)==target_label);
        
        target_nei = neighbor_ins(neighbor_label==target_label);   %% 邻域中与待考察样本类别一致的样本集合
        dis_nei = setdiff(neighbor_ins, target_nei);        %%邻域中其他类别的样本
        dis_nei_label = dataD(dis_nei)';
        D = length(gamma_ds(1,:));%决策数
        ms_local = zeros(1,D);
        ms_local_temp = zeros(1,D);
        ms_all = zeros(1,D);
        %Bel = zeros(N,D);
 
        %% 计算类别标记与待考察样本一致的邻域对象提供的证据信息
        if isempty(target_nei)
            ms_local(1,target_label_ID) = 0;
            ms_all(1,target_label_ID) = 0;
        else
            ms_local_temp(1,target_label_ID) = 1;
            ms_all(1,target_label_ID) = 1;
            for t = 1:length(target_nei)
                ds_value = w1*exp(-target_gamma_ds*(distance(target_nei(t),i).^2))+w2*(1/(1+exp(-value*Sim(target_nei(t),i))));
                ms_local_temp(1,target_label_ID) =  ms_local_temp(1,target_label_ID)*(1-ds_value);
                ms_all(1,target_label_ID) = ms_all(1,target_label_ID)*(1-ds_value);
            end
            ms_local(1,target_label_ID) = 1-ms_local_temp(1,target_label_ID);
        end
        if ~isempty(dis_nei)   %%为真表示邻域中有不同类别样本  
            dis_diff = unique(dis_nei_label);
           %% 根据每个异质类别分别计算局部、全局证据
            for t = 1:length(dis_diff)
                dis_label_ID = find(gamma_ds(1,:)==dis_diff(t));
                dis_gamma_ds = gamma_ds(2,gamma_ds(1,:)==dis_diff(t));
                ms_local_temp(1,dis_label_ID) = 1;
                ms_all(1,dis_label_ID) = 1;
                for j = 1:length(dis_nei)
                    if dis_nei_label(j) == dis_diff(t)
                        ds_value = w1*exp(-dis_gamma_ds*(distance(dis_nei(j),i).^2))+w2*(1/(1+exp(-value*Sim(dis_nei(j),i))));
                        ms_local_temp(1,dis_label_ID) =  ms_local_temp(1,dis_label_ID)*(1-ds_value);
                        ms_all(1,dis_label_ID) = ms_all(1,dis_label_ID)*(1-ds_value);
                    end
                end
                ms_local(1,dis_label_ID) = 1-ms_local_temp(1,dis_label_ID);
            end
            %% 求D-S证据理论中的K值
            nei_diff = unique(neighbor_label);
            if length(nei_diff) == 1
                Bel(i,dis_label_ID) = ms_local(1,dis_label_ID);
            else
                lost_label = setdiff(gamma_ds(1,:), nei_diff);    %% 邻域中未涉及到的标记
                temp_ms = 0;
                if ~isempty(lost_label)
                    lost_label_index = [];
                    all = ms_all;
                    for ll = 1:length(lost_label)
                        lost_label_index = [lost_label_index, find(gamma_ds(1,:)==lost_label(ll))];
                    end
                    all(lost_label_index) = [];
                    for t= 1:length(nei_diff)
                        temp_all = all;
                        temp_all(t) = [];
                        temp_ms = temp_ms+prod(temp_all);
                    end
                    K = temp_ms+(1-length(nei_diff))*prod(all);
                    for t= 1:length(gamma_ds(1,:))
                        temp_all = ms_all;
                        if gamma_ds(1,t) ~= lost_label
                            temp_all([lost_label_index,t]) = [];
                            Bel(i,t) = (ms_local(t)*prod(temp_all))/K;
                            %Pl(1,t) = Bel(1,t)+ms_all(t);
                        end
                    end
                else   %% 邻域中的类别数和原数据一致
                    for t = 1:length(nei_diff)
                        temp_all = ms_all;
                        temp_all(t) = [];
                        temp_ms = temp_ms+prod(temp_all);
                    end
                    K = temp_ms+(1-length(nei_diff))*prod(ms_all);
                    for t = 1:length(gamma_ds(1,:))
                        temp_all = ms_all;
                        temp_all(t) = [];
                        Bel(i,t) = (ms_local(t)*prod(temp_all))/K;
                        %Pl(1,t) = Bel(1,t)+ms_all(t);
                    end
                end
            end
        else  %% 邻域中所有样本的标记与待考察样本的标记一致
            Bel(i,target_label_ID) = ms_local(1,target_label_ID);
            %Pl(1,target_label_ID) = Bel(1,target_label_ID)+ms_all(1,target_label_ID);
        end
        temp_Bel=Bel;
        temp_Bel(i,target_label_ID)=0;
    else
        target_label=dataD(i);
        target_label_ID=find(gamma_ds(1,:)==target_label);
        Bel(i,target_label_ID)=0;
    end
end
end