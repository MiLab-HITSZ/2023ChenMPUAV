function MPAIMA(Global)
% <multi> <real/binary/permutation>
% Nondominated neighbor immune algorithm
% nA ---  20 --- Size of active population
% nC --- 100 --- Size of clone population

%------------------------------- Reference --------------------------------
% M. Gong, L. Jiao, H. Du, and L. Bo, Multiobjective immune algorithm with
% nondominated neighbor-based selection, Evolutionary Computation, 2008,
% 16(2): 225-255.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
            %% Parameter setting
            nA = 50;
            nC = Global.N;
            nCin = 50;
            %% Generate random population
            B = Global.Initialization();               % Antibody population
            [D,FrontNo_low,FrontNo,rp] = UpdateDominantPopulation(B,Global.N,Global.DM,[ones(1,numel(B))]);	% Dominant population
            MaxGen = Global.evaluation/Global.evaluated;
            nowGen = 1;
            stageGen = 0;
            rplist=[];
            rplist(nowGen)=rp;
            p1(nowGen) = 0.9*1/(1+exp(20*(nowGen/MaxGen-0.3)));
            rpmean = mean(rplist(max(end-4,1):end));
            if rpmean<0.1
                stageGen = stageGen+1;
            else
                stageGen = 0;
            end
            p3(nowGen) = 0.9*(1/(1+exp(-20*(stageGen/MaxGen-0.25))));
            %% Optimization
            while Global.NotTermination(D)
                Activeindex = 1:min(min(nA,length(D)));
                Cloneindex = 1:min(min(nCin,length(D)));

                A  = D(Activeindex);            % Active population
                C  = Cloning(D(Cloneindex),nC,FrontNo_low(Cloneindex,:),FrontNo(Cloneindex));                     % Clone population
                C1 = C;
                for n=1:length(C)
                    randnum = rand;
                        if randnum<p1(nowGen) 
                            R = randperm(length(A),4);
                            while(sum(abs(A(R(1)).objs-A(R(2)).objs))==0)
                                R = randperm(length(A),4);
                            end
                            C1(n) = OperatorDE2(C(n),A(R(1)),A(R(2)),A(R(3)),A(R(4)),{0.9,0.7,1,20});
                        else
                            R = randperm(length(A),2);
                            while(sum(abs(A(R(1)).objs-A(R(2)).objs))==0)
                                R = randperm(length(A),2);
                            end
                            if randnum>=p1(nowGen) && p1(nowGen)+p3(nowGen) > randnum 
                                C1(n) = OperatorDE(C(n),A(R(1)),A(R(2)),{0.1,0.5,1,20});
                            else
                                C1(n) = OperatorDE(C(n),A(R(1)),A(R(2)),{0.5,0.5,1,20});
                            end
                        end
                end
                
                [D,FrontNo_low,FrontNo,rp]  = UpdateDominantPopulation([D,C1],Global.N,Global.DM,[ones(1,numel(D)) zeros(1,numel(C1))]);
                nowGen = nowGen+1;
                rplist(nowGen)=rp;
                rpmean = mean(rplist(max(end-4,1):end));
                if rpmean<0.1
                    stageGen = stageGen+1;
                else
                    stageGen = 0;
                end
                p1(nowGen) = 0.9*1/(1+exp(20*(nowGen/MaxGen-0.3)));
                p3(nowGen) = 0.9*(1/(1+exp(-20*(stageGen/MaxGen-0.25))));
            end
end