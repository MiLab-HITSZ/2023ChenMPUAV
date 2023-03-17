classdef MPUAV1< PROBLEM
% <problem> <MPMOP>
    methods
        %% Initialization
        function obj = MPUAV1(data)
            obj.Global.M = 4;
            if isempty(obj.Global.D)
                obj.Global.D =data.dim;
            end
            obj.Global.data = data;
            obj.Global.lower    = data.lb;
            obj.Global.upper    = data.ub;
            obj.Global.DM=2;
            obj.Global.encoding = 'real';
            obj.Global.fn = @callfit_1;
            obj.Global.funname={'length','distance','fatal','eco'};

        end
        %% Calculate objective values for each party
        function PopObj = CalObj(obj,PopDec)         
            PopObj = zeros(size(PopDec,1),obj.Global.M);
            for i =1:size(PopDec,1)
                PopObj(i,:)=obj.Global.fn(PopDec(i,:),obj.Global.data);
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj)
            t1=1;
            t2=2;
            [~,true_y1,true_y2]=true_PS(obj.Global.D,'MPMOP1',t1,t2);
            P=[true_y1,true_y2];
        end
        
        %% Sample reference points on Pareto optimal set
        function P = PS(obj)
            t1=1;
            t2=2;
            P=true_PS(obj.Global.D,'MPMOP1',t1,t2);
        end
        
    end
end