
clear;
clc;
close all;
rng default
rng(1)
addpath(genpath(pwd));
%%
data_2();
data.alpha_trace = 60/360*(2*pi); % 偏航角
data.beta_trace = 30/360*(2*pi); % 仰俯角
data.map_size=map_size;
data.P_crash = 3.42 * 10e-4; % 失控概率
data.S_hit=0.0188; % m^2 撞击面积
data.R_I = 0.3;  % 阻力系数
data.R_vf = 0.27; % 汽车风险
data.alpha=10^6; % J 致死动能
data.beta = 100; % J 
data.S_c = 0.5 ; % 遮蔽系数
data.g = 9.8 ; % m/s^2
data.m = 1380 ; % g (DJI Phantom4)
data.rou_a = 1225 ; % g/m^3(大气密度)
data.miu = miu; % 楼高分布参数
data.sigma = sigma; % 楼高分布参数
data.v = 10; % 20m/s
S=[1 1];E=[45 45]; % 起点终点 
data.S = S;
data.E = E;
data.IOT_pos=IOT_pos;
data.minh=bulid_xyz;
data.maxh=141;
Bound = E(1)-S(1);
dim = Bound*2;
data.Bound = Bound;
data.map_step=map_step;
data.populations_risk=populations_risk;
data.road_risk=road_risk;
GLOBAL
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
RANDSEED=10 ;
count=0;strcell={'(a)','(b)','(c)','(d)'};
for h = 30:30:120
    count=count+1;
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
    title([strcell{count} '  At an altitude of ' num2str(h) ' m'],'Position',[26,-12]);
end
%%
problemList={@MPUAV1,@MPUAV2,@MPUAV3,@MPUAV4,@MPUAV5,@MPUAV6};
maxiterList={40000,40000,40000,40000,40000,40000};
%problemList={@MPUAV1};
%maxiterList={50000};
data.lb = [ones(1,dim/2-1).*-2 ones(1,dim/2+1).*0];
data.ub = [ones(1,dim/2-1).*ystep ones(1,dim/2+1).*1];
data.dim = dim;
figure(100);
for problemindex=1:numel(problemList)
popnum=105;
maxiter=maxiterList{problemindex};
problem=problemList{problemindex};
newp=problem(data);
funname = newp.Global.funname;
callfit = newp.Global.fn;
%% NSGA2
rng default;rng(RANDSEED);
test_case={@OptAll,problem,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population = MPSELECT(Global.result{2},100,2);
[Population2,FrontNo,~] = NDSELECT(Global.result{2},100);
Population2=Population2(FrontNo==1);
res_nsga.dec=reshape([Population.dec],dim,[])';
res_nsga.obj=reshape([Population.obj],numel(funname),[])';
res_nsga2.dec=reshape([Population2.dec],dim,[])';
res_nsga2.obj=reshape([Population2.obj],numel(funname),[])';
%[f,t]=callfit(Population.dec(1,:),data);
% figure()
% hold on 
% colormap('colorcube')
% bar3(bulid_xyz);
% plot(t(:,1),t(:,2),t(:,3),'black-.','LineWidth',1.5);
% xlim([0 100])
% ylim([0 100])
% zlim([0.1 max(bulid_height)+10])
% xlabel('X(50m)');
% ylabel('Y(50m)');
% zlabel('Z(m)');
%%
figure()
title('OptAll')
subplot(1,2,1)
hold on
if numel(funname)==4
plot(res_nsga.obj(:,1),res_nsga.obj(:,2),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
plot(res_nsga2.obj(:,1),res_nsga2.obj(:,2),'MarkerFaceColor','#9FB7A0','Marker','o','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
h=legend({'Common Pareto Solutions','Non-Common Pareto Solutions'});
set(h,'FontName','Times New Roman','FontSize',15)
title('OptALL\_Part1(Efficiency)','FontSize',15)
subplot(1,2,2)
hold on
plot(res_nsga.obj(:,3),res_nsga.obj(:,4),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
plot(res_nsga2.obj(:,3),res_nsga2.obj(:,4),'MarkerFaceColor','#9FB7A0','Marker','o','LineStyle','none');
title('OptALL\_Part2(Safety)','FontSize',15)
xlabel(funname{3},'FontSize',15)
ylabel(funname{4},'FontSize',15)
h=legend({'Common Pareto Solutions','Non-Common Pareto Solutions'});
set(h,'FontName','Times New Roman','FontSize',15)
else
plot3(res_nsga.obj(:,1),res_nsga.obj(:,2),res_nsga.obj(:,3),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
plot3(res_nsga2.obj(:,1),res_nsga2.obj(:,2),res_nsga2.obj(:,3),'MarkerFaceColor','#9FB7A0','Marker','o','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
zlabel(funname{3},'FontSize',15)
h=legend({'Common Pareto Solutions','Non-Common Pareto Solutions'});
set(h,'FontName','Times New Roman','FontSize',15)
title('OptALL\_Part1(Efficiency)','FontSize',15)
subplot(1,2,2)
hold on
plot3(res_nsga.obj(:,4),res_nsga.obj(:,5),res_nsga.obj(:,6),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
plot3(res_nsga2.obj(:,4),res_nsga2.obj(:,5),res_nsga2.obj(:,6),'MarkerFaceColor','#9FB7A0','Marker','o','LineStyle','none');
title('OptALL\_Part2(Safety)','FontSize',15)
xlabel(funname{4},'FontSize',15)
ylabel(funname{5},'FontSize',15)
zlabel(funname{6},'FontSize',15)
h=legend({'Common Pareto Solutions','Non-Common Pareto Solutions'});
set(h,'FontName','Times New Roman','FontSize',15)
end
%%
[f,t]=callfit(res_nsga.dec(1,:),data);
plotres(bulid_xyz,bulid_height,t,'NSGA2 UAV Trace');
plot3(data.S(1),data.S(2),t(find(t(:,1)==data.S(1)),3),'redp','MarkerFaceColor','red','MarkerSize',20);
plot3(data.E(1),data.E(2),t(find(t(:,1)==data.E(1)),3),'redp','MarkerFaceColor','red','MarkerSize',20);
%% MPNDS
rng default;rng(RANDSEED);
test_case={@OptMPNDS,problem,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
res_mpnds.dec=reshape([Population.dec],dim,[])';
res_mpnds.obj=reshape([Population.obj],numel(funname),[])';
%%
figure()
subplot(1,2,1)
if numel(funname)==4
plot(res_mpnds.obj(:,1),res_mpnds.obj(:,2),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
title('OptMPNDS\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot(res_mpnds.obj(:,3),res_mpnds.obj(:,4),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{3},'FontSize',15)
ylabel(funname{4},'FontSize',15)
title('OptMPNDS\_Part2(safety)')
else
plot3(res_mpnds.obj(:,1),res_mpnds.obj(:,2),res_mpnds.obj(:,3),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
zlabel(funname{3},'FontSize',15)
title('OptMPNDS\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot3(res_mpnds.obj(:,4),res_mpnds.obj(:,5),res_mpnds.obj(:,6),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{4},'FontSize',15)
ylabel(funname{5},'FontSize',15)
zlabel(funname{6},'FontSize',15)
title('OptMPNDS\_Part2(safety)')
end
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
%% MPNDS2
rng default;rng(RANDSEED);
test_case={@OptMPNDS2,problem,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
res_mpnds2.dec=reshape([Population.dec],dim,[])';
res_mpnds2.obj=reshape([Population.obj],numel(funname),[])';

%%
figure()
subplot(1,2,1)
if numel(funname)==4
plot(res_mpnds2.obj(:,1),res_mpnds2.obj(:,2),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
title('OptMPNDS2\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot(res_mpnds2.obj(:,3),res_mpnds2.obj(:,4),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{3},'FontSize',15)
ylabel(funname{4},'FontSize',15)
title('OptMPNDS2\_Part2(safety)')
else
plot3(res_mpnds2.obj(:,1),res_mpnds2.obj(:,2),res_mpnds2.obj(:,3),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
zlabel(funname{3},'FontSize',15)
title('OptMPNDS2\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot3(res_mpnds2.obj(:,4),res_mpnds2.obj(:,5),res_mpnds2.obj(:,6),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{4},'FontSize',15)
ylabel(funname{5},'FontSize',15)
zlabel(funname{6},'FontSize',15)
title('OptMPNDS2\_Part2(safety)')
end
legend({'Multiparty Pareto Optimality'},'Location','north')
%% MPNNIA
rng default;rng(RANDSEED);
test_case={@MPNNIA,problem,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
MPNNIA.dec=reshape([Population.dec],dim,[])';
MPNNIA.obj=reshape([Population.obj],numel(funname),[])';
figure()
subplot(1,2,1)
if numel(funname)==4
plot(MPNNIA.obj(:,1),MPNNIA.obj(:,2),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
title('MPNNIA\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot(MPNNIA.obj(:,3),MPNNIA.obj(:,4),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{3},'FontSize',15)
ylabel(funname{4},'FontSize',15)
title('MPNNIA\_Part2(safety)')
else
plot3(MPNNIA.obj(:,1),MPNNIA.obj(:,2),MPNNIA.obj(:,3),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
zlabel(funname{3},'FontSize',15)
title('MPNNIA\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot3(MPNNIA.obj(:,4),MPNNIA.obj(:,5),MPNNIA.obj(:,6),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{4},'FontSize',15)
ylabel(funname{5},'FontSize',15)
zlabel(funname{6},'FontSize',15)
title('OptMPNDS\_Part2(safety)')
end
legend({'Multiparty Pareto Optimality'},'Location','north')
%% MPHEIA
rng default;rng(RANDSEED);
%popnum = 210;
%maxiter=20000;
test_case={@MPHEIA,problem,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
MPHEIA.dec=reshape([Population.dec],dim,[])';
MPHEIA.obj=reshape([Population.obj],numel(funname),[])';
figure()
subplot(1,2,1)
if numel(funname)==4
plot(MPHEIA.obj(:,1),MPHEIA.obj(:,2),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
title('MPHEIA\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot(MPHEIA.obj(:,3),MPHEIA.obj(:,4),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{3},'FontSize',15)
ylabel(funname{4},'FontSize',15)
title('MPHEIA\_Part2(safety)')
else
plot3(MPHEIA.obj(:,1),MPHEIA.obj(:,2),MPHEIA.obj(:,3),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
zlabel(funname{3},'FontSize',15)
title('MPHEIA\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot3(MPHEIA.obj(:,4),MPHEIA.obj(:,5),MPHEIA.obj(:,6),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{4},'FontSize',15)
ylabel(funname{5},'FontSize',15)
zlabel(funname{6},'FontSize',15)
title('MPHEIA\_Part2(safety)')
end
legend({'Multiparty Pareto Optimality'},'Location','north')

%% MPAIMA
rng default;rng(RANDSEED);
%popnum = 210;
%maxiter=20000;
test_case={@MPAIMA,problem,popnum,1,1,maxiter,dim};
for i =1:numel(test_case)/7
var={'-algorithm',test_case{i,1},'-problem',test_case{i,2},'-N',test_case{i,3},'-save',test_case{i,4},'-run',test_case{i,5}, ...,
                        '-evaluation',test_case{i,6},'-D',test_case{i,7},'-data',data};
Global = GLOBAL(var{:});
Global.Start();
end
Population= MPSELECT(Global.result{2},100,2);
MPAIMA.dec=reshape([Population.dec],dim,[])';
MPAIMA.obj=reshape([Population.obj],numel(funname),[])';
figure()
subplot(1,2,1)
if numel(funname)==4
plot(MPAIMA.obj(:,1),MPAIMA.obj(:,2),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
title('MPAIMA\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot(MPAIMA.obj(:,3),MPAIMA.obj(:,4),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{3},'FontSize',15)
ylabel(funname{4},'FontSize',15)
title('MPAIMA\_Part2(safety)')
else
plot3(MPAIMA.obj(:,1),MPAIMA.obj(:,2),MPAIMA.obj(:,3),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
zlabel(funname{3},'FontSize',15)
title('MPAIMA\_Part1(efficiency)')
h=legend({'Multiparty Pareto Optimality'},'Location','north');
set(h,'FontName','Times New Roman','FontSize',15)
subplot(1,2,2)
plot3(MPAIMA.obj(:,4),MPAIMA.obj(:,5),MPAIMA.obj(:,6),'MarkerFaceColor','#874526','Marker','s','LineStyle','none');
xlabel(funname{4},'FontSize',15)
ylabel(funname{5},'FontSize',15)
zlabel(funname{6},'FontSize',15)
title('MPAIMA\_Part2(safety)')
end
legend({'Multiparty Pareto Optimality'},'Location','north')
%%
[~,mi] = min(MPAIMA.obj(:,1));
[f,t]=callfit(MPAIMA.dec(mi,:),data);
plotres(bulid_xyz,bulid_height,t,'MPHEIA UAV Trace')
plot3(data.S(1),data.S(2),t(find(t(:,1)==data.S(1)),3),'redp','MarkerFaceColor','red','MarkerSize',20);
plot3(data.E(1),data.E(2),t(find(t(:,1)==data.E(1)),3),'redp','MarkerFaceColor','red','MarkerSize',20);
plot3(IOT_pos(:,1),IOT_pos(:,2),IOT_pos(:,3).*100,'blueo','MarkerFaceColor','blue','MarkerSize',20)
plot3(t(IOT_pos(:,1),1),t(IOT_pos(:,1),2),t(IOT_pos(:,1),3),'reds','MarkerFaceColor','red','MarkerSize',20)
%%
figure(100)
subplot(2,6,(problemindex-1)*2+1)
hold on
if numel(funname)==4
plot(res_nsga.obj(:,1),res_nsga.obj(:,2),'MarkerFaceColor','#9FB7A0','Marker','*','LineStyle','none');
plot(res_mpnds.obj(:,1),res_mpnds.obj(:,2),'MarkerFaceColor','#26A69A','Marker','^','LineStyle','none');
plot(res_mpnds2.obj(:,1),res_mpnds2.obj(:,2),'MarkerFaceColor','#f7e833','Marker','o','LineStyle','none');
plot(MPNNIA.obj(:,1),MPNNIA.obj(:,2),'MarkerFaceColor','black','Marker','p','LineStyle','none');
plot(MPHEIA.obj(:,1),MPHEIA.obj(:,2),'MarkerFaceColor','#7C2587','Marker','v','LineStyle','none');
plot(MPAIMA.obj(:,1),MPAIMA.obj(:,2),'MarkerFaceColor','blue','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
title('Efficiency DM','FontSize',15)
legend({'NSGA2','MPNDS','MPNDS2','BPNNIA','BPHEIA','BPAIMA'},'Location','best','FontSize',15);
subplot(2,6,(problemindex-1)*2+2)
hold on
plot(res_nsga.obj(:,3),res_nsga.obj(:,4),'MarkerFaceColor','#9FB7A0','Marker','*','LineStyle','none');
plot(res_mpnds.obj(:,3),res_mpnds.obj(:,4),'MarkerFaceColor','#26A69A','Marker','^','LineStyle','none');
plot(res_mpnds2.obj(:,3),res_mpnds2.obj(:,4),'MarkerFaceColor','#f7e833','Marker','o','LineStyle','none');
plot(MPNNIA.obj(:,3),MPNNIA.obj(:,4),'MarkerFaceColor','black','Marker','p','LineStyle','none');
plot(MPHEIA.obj(:,3),MPHEIA.obj(:,4),'MarkerFaceColor','#7C2587','Marker','v','LineStyle','none');
plot(MPAIMA.obj(:,3),MPAIMA.obj(:,4),'MarkerFaceColor','blue','Marker','s','LineStyle','none');
legend({'NSGA2','MPNDS','MPNDS2','BPNNIA','BPHEIA','BPAIMA'},'Location','best','FontSize',15);
xlabel(funname{3},'FontSize',15)
ylabel(funname{4},'FontSize',15)
title('Safety DM','FontSize',15)
else
plot3(res_nsga.obj(:,1),res_nsga.obj(:,2),res_nsga.obj(:,3),'MarkerFaceColor','#9FB7A0','Marker','*','LineStyle','none');
plot3(res_mpnds.obj(:,1),res_mpnds.obj(:,2),res_mpnds.obj(:,3),'MarkerFaceColor','#26A69A','Marker','^','LineStyle','none');
plot3(res_mpnds2.obj(:,1),res_mpnds2.obj(:,2),res_mpnds2.obj(:,3),'MarkerFaceColor','#f7e833','Marker','o','LineStyle','none');
plot3(MPNNIA.obj(:,1),MPNNIA.obj(:,2),MPNNIA.obj(:,3),'MarkerFaceColor','black','Marker','p','LineStyle','none');
plot3(MPHEIA.obj(:,1),MPHEIA.obj(:,2),MPHEIA.obj(:,3),'MarkerFaceColor','#7C2587','Marker','s','LineStyle','none');
plot3(MPAIMA.obj(:,1),MPAIMA.obj(:,2),MPAIMA.obj(:,3),'MarkerFaceColor','#7C2587','Marker','s','LineStyle','none');
xlabel(funname{1},'FontSize',15)
ylabel(funname{2},'FontSize',15)
zlabel(funname{3},'FontSize',15)
title('Efficiency DM','FontSize',15)
legend({'NSGA2','MPNDS','MPNDS2','BPNNIA','BPHEIA','BPAIMA'},'Location','best','FontSize',15);
subplot(2,6,(problemindex-1)*2+2)
hold on
plot3(res_nsga.obj(:,4),res_nsga.obj(:,5),res_nsga.obj(:,6),'MarkerFaceColor','#9FB7A0','Marker','*','LineStyle','none');
plot3(res_mpnds.obj(:,4),res_mpnds.obj(:,5),res_mpnds.obj(:,6),'MarkerFaceColor','#26A69A','Marker','^','LineStyle','none');
plot3(res_mpnds2.obj(:,4),res_mpnds2.obj(:,5),res_mpnds2.obj(:,6),'MarkerFaceColor','#f7e833','Marker','o','LineStyle','none');
plot3(MPNNIA.obj(:,4),MPNNIA.obj(:,5),MPNNIA.obj(:,6),'MarkerFaceColor','black','Marker','p','LineStyle','none');
plot3(MPHEIA.obj(:,4),MPHEIA.obj(:,5),MPHEIA.obj(:,6),'MarkerFaceColor','#7C2587','Marker','s','LineStyle','none');
plot3(MPAIMA.obj(:,1),MPAIMA.obj(:,2),MPAIMA.obj(:,3),'MarkerFaceColor','#7C2587','Marker','s','LineStyle','none');
legend({'NSGA2','MPNDS','MPNDS2','BPNNIA','BPHEIA','BPAIMA'},'Location','best','FontSize',15);
xlabel(funname{4},'FontSize',15)
ylabel(funname{5},'FontSize',15)
zlabel(funname{6},'FontSize',15)
title('Safety DM','FontSize',15)
end
%%
allsol = [res_nsga.obj;res_mpnds.obj;res_mpnds2.obj;MPNNIA.obj;MPHEIA.obj];
nsga_hv = MPHV(res_nsga.obj,allsol,2)
mpnds_hv = MPHV(res_mpnds.obj,allsol,2)
mpnds2_hv = MPHV(res_mpnds2.obj,allsol,2)
MPNNIA_hv = MPHV(MPNNIA.obj,allsol,2)
MPHEIA_hv = MPHV(MPHEIA.obj,allsol,2)
MPAIMA_hv = MPHV(MPAIMA.obj,allsol,2)
end
%%
function plotres(bulid_xyz,bulid_height,t,titlestr)
figure()
hold on 
colormap('colorcube')
remove_empty_bars(bar3(bulid_xyz,1));
plot3(t(:,1),t(:,2),t(:,3),'Color','#A2AB28','LineStyle','-','Marker','o','LineWidth',1.5);
xlim([0 50])
ylim([0 50])
zlim([0 max(bulid_height)+10])
xlabel('X(unit/100m)');
ylabel('Y(unit/100m)');
zlabel('Z(unit/1m)');
title(titlestr)
end 

