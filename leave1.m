function [K] = leave1(Data,rho)
DataC=Data(:,1:end-1);
[Data_num,~]=size(DataC);
Raw_data=DataC;
%参数调节
opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% termination criterion
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=100;   % maximum number of iterations

% normalization
opts.nFlag=0;       % without normalization

% regularization
opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
%opts.rsL2=0.01;     % the squared two norm term
% rho=0.02;
for i=1:Data_num
    current_sample=DataC(i,:);
    DataC(i,:)=[];
    A=DataC';%转置摆放
    y=current_sample';
    [W1, ~, ~]= LeastR(A, y, rho, opts);
%     m=sum((W1~=0));
%     m=sum((W1~=0));
    K(:,i)=W1;
    DataC=Raw_data;
end
end

