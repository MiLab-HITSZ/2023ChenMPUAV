function C = Cloning2(A,nC,FrontNo_low,FrontNo)
% Proportional cloning

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the crowding distance of each solution
    [PopObj,~,Loc] = unique(A.objs,'rows');
    M = size(PopObj,2);
    DM = size(FrontNo_low,2);
    nC = nC/DM;
    for d =1:DM
        CrowdDis = CrowdingDistance(PopObj(:,1+(d-1)*(M/DM):d*(M/DM)));
        CrowdDis = CrowdDis(Loc);
        if all(CrowdDis==inf)
            CrowdDis = ones(size(CrowdDis));
        else
            CrowdDis(CrowdDis==inf) = 2*max(CrowdDis(CrowdDis~=inf));
        end
        %% Clone
        q = ceil(nC*CrowdDis/sum(CrowdDis));
        C = [];
        for i = 1 : length(A)
            C = [C,repmat(A(i),1,q(i))];
        end
    end
end