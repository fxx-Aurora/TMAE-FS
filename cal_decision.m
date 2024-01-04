function [decisionclass, different] = cal_decision(dataD)
%cal_decision 函数返回决策属性的决策类
%Input：dataD 决策属性集
%output: decisionclass 决策类集合，为一n*m的矩阵，n为对象数，m为决策类

different = unique(dataD);
[m,~] = size(different);%决策类个数
[n,~] = size(dataD);%样本个数
decision = zeros(n,m);
for i = 1:m
    for j = 1:n
        if dataD(j) == different(i)
            decision(j,i) = 1;
        end
    end
end
decisionclass = sparse(decision);%去除大量非0样本
%decisionclass = sparse(bsxfun(@eq,dataD(:,1),rot90(dataD(:,1))));
clear m n dataD i j decision tmp tmpdecision
end