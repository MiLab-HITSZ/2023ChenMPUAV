function Fitness = CalD_SPEASDE(PopObjALL,Global)
% Calculate the fitness of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    M = Global.M;
    N = size(PopObjALL,1);
    DM = Global.DM;
    D = zeros(N,DM);

    for dms =1:Global.DM
        PopObj = PopObjALL(:,1+(dms-1)*M/DM:dms*M/DM);
        
               
        %% Calculate the shifted distance between each two solutions
        Distance = inf(N);
        for i = 1 : N
            SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
            for j = [1:i-1,i+1:N]
                Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
            end
        end
        
        %% Calculate D(i)
        Distance = sort(Distance,2);
        D(:,dms) = 1./(Distance(:,floor(sqrt(N)))+2);
    end
    %% Calculate the fitnesses
    Fitness = sum(D');
end