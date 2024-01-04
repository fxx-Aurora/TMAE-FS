function [Approximation_quality]=Depdency_TMVN(data,rh)
[N,~]=size(data);
dataD=data(:,end);
[~,Neighbor,~,~,~]=Neigh(data,rh);

[DC,~]=cal_decision(dataD);
[~,~,~,Cardinal]=cal_neighbor_low(DC,Neighbor); %计算下近似，正域
Approximation_quality=sum(Cardinal)/N; % 依赖度
end