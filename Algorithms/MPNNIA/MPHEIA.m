function MPHEIA(Global)
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
            [D,FrontNo_low,FrontNo] = UpdateDominantPopulation(B,Global.N,Global.DM);	% Dominant population

            %% Optimization
            while Global.NotTermination(D)
                Activeindex = 1:min(min(nA,length(D)));
                Cloneindex = 1:min(min(nCin,length(D)));
                A  = D(Activeindex);            % Active population
                try
                C  = Cloning(D(Cloneindex),nC,FrontNo_low(Cloneindex,:),FrontNo(Cloneindex));                     % Clone population
                catch
                    disp('error');
                end
                C1 = C;
                for n=1:length(C)
                    if rand<0
                        R = randperm(length(A),2);
                        while(sum(abs(A(R(1)).objs-A(R(2)).objs))==0)
                            R = randperm(length(A),2);
                        end
                        C1(n) = OperatorDE(C(n),A(R(1)),A(R(2)));
                    else
                        R = randi(length(A));
                        while(sum(abs(C(n).objs-A(R).objs))==0)
                            R = randi(length(A));
                        end
                        C1(n) = OperatorGAhalf([C(n),A(R)]);
                    end
                end
                
                [D,FrontNo_low,FrontNo]  = UpdateDominantPopulation([D,C1],Global.N,Global.DM);
            end
end