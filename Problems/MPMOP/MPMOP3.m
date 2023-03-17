classdef MPMOP3< PROBLEM
% <problem> <MPMOP>
    methods
        %% Initialization
        function obj = MPMOP3()
            obj.Global.M = 4;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    =ones(1,obj.Global.D).*(-1);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.lower(1)=0;
            obj.Global.DM=2;
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values for each party
        function PopObj = CalObj(obj,PopDec)
            M=obj.Global.M;
            PopObj(:,[1:M/2])=MPMOP_Value('MPMOP3', PopDec(), 0);
            PopObj(:,[M/2+1:M])=MPMOP_Value('MPMOP3', PopDec(),pi/2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj)
            t1=0;
            t2=pi/2;
            [~,true_y1,true_y2]=true_PS(obj.Global.D,'MPMOP3',t1,t2);
            P=[true_y1,true_y2];
        end
        
        %% Sample reference points on Pareto optimal set
        function P = PS(obj)
            t1=0;
            t2=pi/2;
            P=true_PS(obj.Global.D,'MPMOP3',t1,t2);
        end
     
        
    end
end