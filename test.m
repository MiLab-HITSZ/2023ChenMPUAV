clear;
clc;
close all;
rng default
rng(11)
addpath(genpath(pwd));
%%
data_2();
data.alpha_trace = 60/360*(2*pi); % 偏航角
data.beta_trace = 45/360*(2*pi); % 仰俯角
data.map_size=map_size;
data.P_crash = 3.42 * 10e-4; % 失控概率
data.S_hit=0.0188; % m^2 撞击面积
data.R_I = 0.3;  % 阻力系数
data.R_vf = 0.27; % 汽车风险
data.alpha=10^6; % J 致死动能
data.beta = 100; % J 
data.S_c = 0.5 ; % 遮蔽系数
data.g = 9.8 ; % m/s^2
data.IOT_pos=IOT_pos;
data.m = 1380 ; % g (DJI Phantom4)
data.rou_a = 1225 ; % g/m^3(大气密度)
data.miu = miu; % 楼高分布参数
data.sigma = sigma; % 楼高分布参数
data.v = 20; % 20m/s
S=[1 1];E=[49 49]; % 起点终点 
data.S = S;
data.E = E;
data.minh=bulid_xyz;
data.maxh=141;
Bound = E(1)-S(1);
dim = Bound*2;
data.Bound = Bound;
data.map_step=map_step;
data.populations_risk=populations_risk;
data.road_risk=road_risk;
%% pre-cal
ystep = 3;
pbase = ystep+1;
for i = 1:2*ystep+1
    pi = i - pbase;
    can=[];
    for j = -ystep:1:ystep
        if acos([1,pi]*[1,j]'/sqrt(1+pi^2)/sqrt(1+(j)^2))<=data.alpha_trace
            can=[can j];
        end
    end
    canselect{i}=can;
end
data.canselect = canselect;
data.canselectp = pbase;
%%
tiledlayout(2,2);
for h = 30:30:120
    nexttile;
    Risk_map = zeros(map_size);
    Riskproperty_map = zeros(map_size);
    for i=1:map_size(1)
        for j =1:map_size(2)
            Risk_map(i,j)=Risk_map(i,j)+getC_Risk(getR_pf(getV(h,data),data),populations_risk(i,j),data);
            Risk_map(i,j)=Risk_map(i,j)+getC_Risk(data.R_vf,road_risk(i,j),data);
        end
    end
    colormap('jet')
    contourf(Risk_map)
    colorbar;
    title(['h=' num2str(h) 'm,' ' Risk of property=' num2str(getC_rpd(h,data))]);
end
%%
problemList={@MPUAV1,@MPUAV2,@MPUAV3,@MPUAV4,@MPUAV5,@MPUAV6};
maxiterList={40000,40000,40000,40000,40000,40000};
problemMean=zeros(numel(problemList),6);
problemStd=zeros(numel(problemList),6);
data.lb = [ones(1,dim/2-1).*-1 ones(1,dim/2+1).*0];
data.ub = [ones(1,dim/2-1).*ystep ones(1,dim/2+1).*1];
data.dim = dim;
temp.dec=0;
temp.obj=0;
for problemIndex= 1:numel(problemList)
    TT=30;
score=zeros(TT,6);
% res_nsga=repmat(temp,1,numel(problemList));
% res_nsga2=repmat(temp,1,numel(problemList));
% res_mpnds=repmat(temp,1,numel(problemList));
% res_mpnds2=repmat(temp,1,numel(problemList));
% MPNNIA=repmat(temp,1,numel(problemList));
% MPHEIA=repmat(temp,1,numel(problemList));
% MPAIMA=repmat(temp,1,numel(problemList));
%RS=1:31;RS(30)=[];
parfor testtimes = 10:TT
close all;
RANDSEED=testtimes;
testfit = problemList{problemIndex};
popnum=105;
maxiter=maxiterList{problemIndex};
%% NSGA2
rng default;rng(RANDSEED);
test_case={@OptAll,testfit,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population = MPSELECT(Global.result{2},100,2);
[Population2,FrontNo,~] = NDSELECT(Global.result{2},100);
Population2=Population2(FrontNo==1);
res_nsgadec=reshape([Population.dec],dim,[])';
res_nsgaobj=reshape([Population.obj],size(Population(1).obj,2),[])';
res_nsga2dec=reshape([Population2.dec],dim,[])';
res_nsga2obj=reshape([Population2.obj],size(Population(1).obj,2),[])';
%% MPNDS
rng default;rng(RANDSEED);
test_case={@OptMPNDS,testfit,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
res_mpndsdec=reshape([Population.dec],dim,[])';
res_mpndsobj=reshape([Population.obj],size(Population(1).obj,2),[])';
%% MPNDS2 
rng default;rng(RANDSEED);
test_case={@OptMPNDS2,testfit,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
res_mpnds2dec=reshape([Population.dec],dim,[])';
res_mpnds2obj=reshape([Population.obj],size(Population(1).obj,2),[])';
%% MPNNIA
rng default;rng(RANDSEED);
%popnum = 210;
%maxiter=20000;
test_case={@MPNNIA,testfit,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
MPNNIAdec=reshape([Population.dec],dim,[])';
MPNNIAobj=reshape([Population.obj],size(Population(1).obj,2),[])';
%% MPHEIA
rng default;rng(RANDSEED);
%popnum = 210;
%maxiter=20000;
test_case={@MPHEIA,testfit,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
MPHEIAdec=reshape([Population.dec],dim,[])';
MPHEIAobj=reshape([Population.obj],size(Population(1).obj,2),[])';
%% MPAIMA
rng default;rng(RANDSEED);
%popnum = 210;
%maxiter=20000;
test_case={@MPAIMA,testfit,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
MPAIMAdec=reshape([Population.dec],dim,[])';
MPAIMAobj=reshape([Population.obj],size(Population(1).obj,2),[])';
%%
allsol = [res_nsgaobj;res_mpndsobj;res_mpnds2obj;MPNNIAobj;MPHEIAobj;MPAIMAobj];
nsga_hv = MPHV(res_nsgaobj,allsol,2);
mpnds_hv = MPHV(res_mpndsobj,allsol,2);
mpnds2_hv = MPHV(res_mpnds2obj,allsol,2);
MPNNIA_hv = MPHV(MPNNIAobj,allsol,2);
MPHEIA_hv = MPHV(MPHEIAobj,allsol,2);
MPAIMA_hv = MPHV(MPAIMAobj,allsol,2);
score(testtimes,:)=[nsga_hv mpnds_hv mpnds2_hv MPNNIA_hv MPHEIA_hv MPAIMA_hv];
end
problemMean(problemIndex,:)= mean(score)
clear var
problemStd(problemIndex,:)=var(score)
end
% problemMean =
% 
%     0.0483    0.0682    0.0696    0.0519    0.0673    0.0807
%     0.0533    0.1022    0.1051    0.0726    0.1102    0.1106
%     0.0519    0.0843    0.0785    0.0610    0.0899    0.0914
%     0.0494    0.0589    0.0523    0.0424    0.0539    0.0663
%     0.0439    0.0744    0.0771    0.0518    0.0775    0.0788
%     0.0369    0.0615    0.0572    0.0403    0.0598    0.0640
% 
% 
% problemStd =
% 
%     0.0003    0.0002    0.0002    0.0003    0.0004    0.0004
%     0.0009    0.0007    0.0006    0.0007    0.0009    0.0010
%     0.0006    0.0004    0.0004    0.0005    0.0005    0.0006
%     0.0007    0.0003    0.0004    0.0005    0.0005    0.0004
%     0.0006    0.0007    0.0007    0.0009    0.0012    0.0008
%     0.0004    0.0004    0.0003    0.0003    0.0003    0.0005
% problemMean =
% 
%     0.0485    0.0681    0.0689    0.0669    0.0698    0.0800
%     0.0530    0.1021    0.1038    0.1042    0.1013    0.1108
%     0.0505    0.0817    0.0757    0.0830    0.0795    0.0891
%     0.0488    0.0580    0.0513    0.0445    0.0448    0.0654
%     0.0451    0.0770    0.0793    0.0636    0.0657    0.0803
%     0.0362    0.0610    0.0565    0.0539    0.0533    0.0641
% 
% 
% problemStd =
% 
%     0.0003    0.0002    0.0003    0.0003    0.0003    0.0004
%     0.0009    0.0007    0.0006    0.0012    0.0008    0.0009
%     0.0005    0.0004    0.0004    0.0004    0.0004    0.0006
%     0.0007    0.0003    0.0004    0.0004    0.0005    0.0004
%     0.0007    0.0006    0.0007    0.0009    0.0010    0.0008
%     0.0004    0.0005    0.0004    0.0004    0.0004    0.0005