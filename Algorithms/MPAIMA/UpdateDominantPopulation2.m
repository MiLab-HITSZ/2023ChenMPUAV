function [D,rp] = UpdateDominantPopulation2(D,N,DM,FLAG,Global)
% Update the dominant population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    objs = D.objs;
    [~,uindex,~] = unique(objs,'rows');
    D = D(uindex);

    Fitness = CalFitness_SPEASDE(D.objs,Global);
    FLAG = FLAG(uindex);
    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(D(Next).objs,sum(Next)-N,DM,Global.M);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    D = D(Next);
    F = Fitness(Next);
    rp = 1-sum(FLAG(Next))/sum(FLAG);
end

function Del = Truncation(PopObjALL,K,DM,M)
% Select part of the solutions by truncation

    %% Truncation
    DistanceALL=cell(DM);
    Distance = 0;
    for dms =1:DM
        PopObj = PopObjALL(:,1+(dms-1)*M/DM:dms*M/DM);
        DistanceALL{dms} = pdist2(PopObj,PopObj);
        DistanceALL{dms}(logical(eye(length(DistanceALL{dms})))) = inf;
        Distance = Distance + DistanceALL{dms};
    end
    
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end