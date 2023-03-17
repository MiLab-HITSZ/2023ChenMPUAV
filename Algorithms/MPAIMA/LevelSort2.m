function [FrontNo,MaxFNo,FrontNo_low] = LevelSort2(varargin)
%NDSort - Do non-dominated sorting by efficient non-dominated sort.
%
%   FrontNo = NDSort(F,s) does non-dominated sorting on F, where F is the
%   matrix of objective values of a set of individuals, and s is the number
%   of individuals to be sorted at least. FrontNo(i) denotes the front
%   number of the i-th individual. The individuals have not been sorted are
%   assigned a front number of inf.
%
%   FrontNo = NDSort(F,C,s) does non-dominated sorting based on constrained
%   domination, where C is the matrix of constraint values of the
%   individuals. In this case, feasible solutions always dominate
%   infeasible solutions, and one infeasible solution dominates another
%   infeasible solution if the former has a smaller overall constraint
%   violation than the latter.
%
%   In particular, s = 1 indicates finding only the first non-dominated
%   front, s = size(F,1)/2 indicates sorting only half the population
%   (which is often used in the algorithm), and s = inf indicates sorting
%   the whole population.
%
%   [FrontNo,K] = NDSort(...) also returns the maximum front number besides
%   inf.
%
%   Example:
%       [FrontNo,MaxFNo] = NDSort(PopObj,1)
%       [FrontNo,MaxFNo] = NDSort(PopObj,PopCon,inf)

%------------------------------- Reference --------------------------------
% [1] X. Zhang, Y. Tian, R. Cheng, and Y. Jin, An efficient approach to
% nondominated sorting for evolutionary multiobjective optimization, IEEE
% Transactions on Evolutionary Computation, 2015, 19(2): 201-213.
% [2] X. Zhang, Y. Tian, R. Cheng, and Y. Jin, A decision variable
% clustering based evolutionary algorithm for large-scale many-objective
% optimization, IEEE Transactions on Evolutionary Computation, 2018, 22(1):
% 97-112.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = varargin{1};
    [N,M]  = size(PopObj);
    if nargin == 2
        nSort  = varargin{2};
    else
        PopCon = varargin{2};
        nSort  = varargin{3};
        Infeasible           = any(PopCon>0,2);
        PopObj(Infeasible,:) = repmat(max(PopObj,[],1),sum(Infeasible),1) + repmat(sum(max(0,PopCon(Infeasible,:)),2),1,M);
    end
    Tabletemp = cell(1,M);
    cumTable = cell(1,M);
    for i=1:M
        Tabletemp{1,i} = cumsum(hist(PopObj(:,i),max(PopObj(:,i))));
        cumTable{1,i} = [0 Tabletemp{1,i}];
    end
    [PopObj,~,Loc] = unique(PopObj,'rows');
    Table          = hist(Loc,1:max(Loc));
    [N,M]          = size(PopObj);
    FrontNo = inf(1,N); 
    MaxFNo         = 0;
    
    RANK_SCORE1 = zeros(size(PopObj));
    RANK_SCORE2 = zeros(size(PopObj));
    for i=1:M
        RANK_SCORE(:,i) = cumTable{i}(PopObj(:,i));
        RANK_SCORE2(:,i) = Tabletemp{i}(PopObj(:,i));
    end
    [~,minindex]=min(sum(RANK_SCORE,2));
    FrontNo(minindex)=1;
    while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
        MaxFNo = MaxFNo + 1;
        for i = 1 : N
            if FrontNo(i) == inf
                Dominated = false;
                for j = 1:N
                    if FrontNo(j) == MaxFNo
                        Dominated=LevelDominated(PopObj(j,:),PopObj(i,:),RANK_SCORE(j,:),...
                            RANK_SCORE(i,:),RANK_SCORE2(j,:),RANK_SCORE2(i,:));
                        if Dominated
                            break;
                        end
                    end
                end
                if ~Dominated
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end
    FrontNo = FrontNo(:,Loc);
    MaxFNo = max(MaxFNo,1);
end
function Dominated=LevelDominated(x,y,xs1,ys1,xs2,ys2)
Dominated = 0;
if any(x<y) && all(x<=y)
    Dominated=1;
    return
end
sxs1=sort(xs1);sys1=sort(ys1);
if sys1(1)~=1 && sxs1(1)~=1 && any(sxs1<sys1) && all(sxs1<sys1)
        Dominated=1;
    return
end
% sxs2=sort(xs2);sys2=sort(ys2);
% if any(sxs2<sys2) && all(sxs2<sys2)
%         Dominated=1;
%     return
% end
end