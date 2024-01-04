function [NDER_value1,ave]=NDER(data,Neighbor)
dataD=data(:,end);
[DC, DC_diff] = cal_decision(dataD);
[min_label,min_index]=min(DC_diff);
if min_label==0
    dataD(find(dataD==0))=max(DC_diff)+1;
    DC_diff(min_index)=max(DC_diff)+1;
end
[N,~]=size(data);
[~,M]=size(DC);
%邻域决策错误率
NDER_value=0;
for i=1:N
    Neighborhood=find(Neighbor(:,i)==1);  %每个样本的邻域集合
    [a,~]=size(Neighborhood);%邻域的个数
    prob=[];
    for j=1:M
        X=find((DC(:,j))==1);
        Y=intersect(Neighborhood,X);
        [b,~]=size(Y);
        prob_temp=b/a;
        prob(DC_diff(j))=prob_temp;
    end
    max_value=max(prob);
    for j=1:M
        if max_value==prob(j)
            pre_label=j;
        end
    end
    if pre_label==dataD(i)
        loos_value=0;
    else
        loos_value=1;
    end
    NDER_value=NDER_value+loos_value;
end
NDER_value=NDER_value/N;
NDER_value1=1-NDER_value; %1-邻域决策错误率
for i=1:N
    temp=0; %计算邻居中的类别和该样本本身的类别相同的个数
    Neighborhood=find(Neighbor(:,i)==1);  %每个样本的邻域集合
    [a,~]=size(Neighborhood);%邻域的个数
    if a<0.01*N
        Nei(i)=0;
    else
        for j=1:a
            if dataD(i)==dataD(Neighborhood(j))
                temp=temp+1;
            end
        end
        Nei(i)=temp/a;
    end
end
ave=mean(Nei);
end