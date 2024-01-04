function [decisionclass, different] = cal_decision(dataD)
%cal_decision �������ؾ������Եľ�����
%Input��dataD �������Լ�
%output: decisionclass �����༯�ϣ�Ϊһn*m�ľ���nΪ��������mΪ������

different = unique(dataD);
[m,~] = size(different);%���������
[n,~] = size(dataD);%��������
decision = zeros(n,m);
for i = 1:m
    for j = 1:n
        if dataD(j) == different(i)
            decision(j,i) = 1;
        end
    end
end
decisionclass = sparse(decision);%ȥ��������0����
%decisionclass = sparse(bsxfun(@eq,dataD(:,1),rot90(dataD(:,1))));
clear m n dataD i j decision tmp tmpdecision
end