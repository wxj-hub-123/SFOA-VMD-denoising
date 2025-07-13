function [xposbest,xposprocess,fvalbest,Curve]=SFOA_VMD(x,Npop,Max_it,lb,ub)
% SFOA_VMD   Using starfish optimization algorithm to calculate the optimal parameters of VMD.
%
%   [xposbest,xposprocess,fvalbest,Curve]=SFOA_VMD(x,Npop,Max_it,lb,ub) Calculate the
%   optimal parameters for VMD, Enable the decomposed IMF to reflect the true 
%   situation of the data. 
%
%   Input:
%       x - input data
%       Npop - The number of starfish
%       Max_it - Maximum Number Of Iterations
%       lb - Parameters upper limit
%       ub - Parameters lower limit 
%
%   Output:
%       xposbest - optimum parameter
%       xposprocess -  The best result of each iteration
%       fvalbest - The fitness function value under optimal parameter conditions
%       Curve - The fitness function value during each iteration process
%
%   License: MIT
%   Part of the SFOA-VMD denoising framework (version 1.0)
%   Repository: https://github.com/anonymous-researcher/SFOA-VMD-denoising
%
%   Author: Anonymous (for peer review)
%
% -------------------------------------------------------------------------
%% Initialization
GP=0.5;     % parameter 
nD=2;
if size(ub,2) == 1
    lb = lb*ones(1,nD); 
    ub = ub*ones(1,nD);
end

fvalbest = inf;
Fitness=zeros(Npop,1);
Curve = zeros(1,Max_it);
Xpos(:,1) = rand(Npop,1).*(ub(1)-lb(1))+lb(1);   
Xpos(:,2) = rand(Npop,1).*(ub(2)-lb(2))+lb(2);  % 位置X初始化 二维 
for i=1:Npop
    [imf,~,~] = vmd(x,'NumIMF',round(Xpos(i,1)),'PenaltyFactor',round(Xpos(i,2)));
    for j=1:round(Xpos(i,1))
        entro(j)=min_envelope_entropy(imf(:,j));
    end
    Fitness(i)=1.0/abs(max(entro)-min(entro));
end
[fvalbest,order]=min(Fitness);
xposbest=Xpos(order,:);

xposprocess(:,1)=xposbest;
newX = zeros(Npop,nD);

%% Evolution
T = 1;
Curve=zeros(Max_it,1);
while T <= Max_it
    theta = pi/2*T./Max_it;  % theta海星手臂弯曲
    tEO = (Max_it-T)/Max_it*cos(theta);  %Et海星能量
    if rand < GP %  勘探阶段 判断r 
        for i = 1:Npop
            if nD > 5
                % for nD is larger than 5 判断维度大于5 分两个形式计算
                jp1 = randperm(nD,5);
                for j = 1:5
                    pm = (2*rand-1)*pi;
                    if rand < GP
                        newX(i,jp1(j)) = Xpos(i,jp1(j)) + pm*(xposbest(jp1(j))-Xpos(i,jp1(j)))*cos(theta);
                    else
                        newX(i,jp1(j)) = Xpos(i,jp1(j)) - pm*(xposbest(jp1(j))-Xpos(i,jp1(j)))*sin(theta);
                    end
                    if newX(i,jp1(j))>ub(jp1(j)) || newX(i,jp1(j))<lb(jp1(j))
                        newX(i,jp1(j)) = Xpos(i,jp1(j));
                    end
                end
            else
                % for nD is not larger than 5 判断维度小于等于5
                jp2 = ceil(nD*rand); % ceil 向上取整 选取维度jp2
                im = randperm(Npop); % 生成1-50的随机数 选取海星的位置
                rand1 = 2*rand-1;  % A1
                rand2 = 2*rand-1;  % A2
                newX(i,jp2) = tEO*Xpos(i,jp2) + rand1*(Xpos(im(1),jp2)-Xpos(i,jp2))+rand2*(Xpos(im(2),jp2)-Xpos(i,jp2)); % 计算海星的新位置
                if newX(i,jp2)>ub(jp2) || newX(i,jp2)<lb(jp2) %判断位置是否超出边界 超出边界则取原位置
                    newX(i,jp2) = Xpos(i,jp2);
                end  
            end
            newX(i,:) = max(min(newX(i,:),ub),lb);  % boundary check 检测值是否在边界内
        end
    else % 海星的开发利用
        df = randperm(Npop,5); % 生成五个 1-50的随机数 随机选取五个海星
        dm(1,:) = xposbest - Xpos(df(1),:);
        dm(2,:) = xposbest - Xpos(df(2),:);
        dm(3,:) = xposbest - Xpos(df(3),:);
        dm(4,:) = xposbest - Xpos(df(4),:);
        dm(5,:) = xposbest - Xpos(df(5),:);  % five arms of starfish 随机五个海星与最佳位置的距离
        for i = 1:Npop
            r1 = rand; r2 = rand;
            kp = randperm(length(df),2); % 选取两个海星
            newX(i,:) = Xpos(i,:) + r1*dm(kp(1),:) + r2*dm(kp(2),:);   % exploitation 更新海星位置
            if i == Npop
                newX(i,:) = exp(-T*Npop/Max_it).*Xpos(i,:);  % regeneration of starfish 海星再生
            end
            newX(i,:) = max(min(newX(i,:),ub),lb);  % boundary check 检测值是否在边界内
        end
    end
    
    % Fitness evaluation  适应度评价
    for i=1:Npop
        [imf,~,~] = vmd(x,'NumIMF',round(newX(i,1)),'PenaltyFactor',round(newX(i,2)));
        %   newFit=kurt_env(imf); %计算峭度 适应度函数
        %   newFit=min_envelope_entropy(imf); %最小包络熵适应度函数
        min1=min_envelope_entropy(imf(:,end));
        max1=min_envelope_entropy(imf(:,1));
        newFit=1.0/abs(max1-min1);
        if newFit<Fitness(i)
            Fitness(i) = newFit;  % 替换适应度
            Xpos(i,:) = newX(i,:); % 替换位置
            if newFit < fvalbest % 判断适应度是否最佳
                fvalbest = Fitness(i); %替换最佳适应度
                xposbest = Xpos(i,:); %替换最佳位置
            end
        end
    end
   
    Curve(T) = fvalbest; %记录本次迭代的适应度
    T = T+1;
    xposprocess(:,T)=xposbest;
end
end