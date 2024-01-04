function [Reduct_location,Location_index]=Heurstic_TMAEFS(data,rho)
%输入：数据集
%输出：Reduct_location：候选属性的位置，Location_index：被选中的属性值为1，其余为0
%启发式约简算法
%初始化数据
[U_ent]=entory(data,rho);
Contional_entory=U_ent;

dataC=data(:,1:end-1);
dataD=data(:,end);
Selected_data=dataC;

[~,Attribute_num]=size(Selected_data);
Reduct_location=-1*ones(1,Attribute_num);  %计算侯选属性位置
ENT=-1*ones(1,Attribute_num);
Location_index=zeros(1,Attribute_num);  %位置索引，被选中的属性为1，其余为0，初始化全为0
First_turnlist=zeros(1,Attribute_num);  %记录属性列表
for i=1:Attribute_num
    Temp_data=[Selected_data,dataD];
    Temp_data(:,i)=[];  %删除对应的属性
    [ENT_num]=entory(Temp_data,rho);
    First_turnlist(1,i)=ENT_num;  %计算内部属性的条件熵（重要度）
end
%%启发式算法计算下近似约简
[~,first_index]=max(First_turnlist); %选取最大值的及其对应的位置
Temp_data=Selected_data(:,first_index);
Location_index(1,first_index)=1;
Reduct_location(1,1)=first_index;
ENT(1,1)=entory([Temp_data,dataD],rho);
ENT_num=ENT(1,1);

for m=2:Attribute_num
    if ENT_num<U_ent
        U_ent=ENT_num;
    end
    LTMP=[];
    for i=1:Attribute_num
        if ismember(i,Reduct_location(1,:))==0
            Temp_data=[Temp_data,Selected_data(:,i),dataD]; %Temp_data只有条件属性，加入一个属性，再加一个决策属性
            [ENT_region1]=entory(Temp_data,rho);
            LTMP(1,end+1)=i; %记录属性
            LTMP(2,end)=ENT_num-ENT_region1;  %计算sig(A∪{a},d-sig(A,d))
            Temp_data(:,end-1:end)=[];%去掉上次的条件属性
        end
    end
    [~,X_index]=max(LTMP(2,:)); %最大的属性重要度
    first_index=LTMP(1,X_index); %记录属性的位置
    %将该属性的属性列，属性标志等信息存储下来
    Reduct_location(1,m)=first_index;
    Location_index(1,first_index)=1;
    Temp_data=[Temp_data,Selected_data(:,first_index),dataD];%这里的Temp_data只有条件属性，加入一个属性，再加一个决策属性
    [ENT_num]=entory(Temp_data,rho);
    ENT(1,m)=ENT_num;
    if ENT_num<U_ent
        U_ent=ENT_num;
    end
    if ENT_num>U_ent&&ENT(1,m-1)<Contional_entory
        Reduct_location=Reduct_location(1,1:m);
        Location_index(1,first_index)=0;
        ENT=ENT(1,1:m-1);
        break
    end
    Temp_data(:,end)=[];
end
%进行属性删减的工作，对于约简集中的每一个属性，删除判断熵如果不降或者上升，那么就将这个元素删除
[~,Reduct_num]=size(Reduct_location);
Target_ENT=ENT(1,end);
for w=1:Reduct_num
    if Reduct_num<=2
        break
    end
    red=Reduct_location;
    red(:,w)=0;
    red(:,red==0)=[];
    Test_data=[Selected_data(:,red),dataD];
    [ENT_num_test]=entory(Test_data,rho);
    if ENT_num_test<Target_ENT
        Reduct_location(:,w)=0;
        ENT(1,end+1)=ENT_num_test;
        Target_ENT=ENT_num_test;
    end
end
Reduct_location(:,Reduct_location==0)=[];
end