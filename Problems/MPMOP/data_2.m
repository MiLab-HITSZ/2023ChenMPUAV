rng(12)
max_altitude = 120;
m_size = [5000 5000]; % m
map_size = [50 50]; % nums
map_mat = zeros(map_size); % mapMat
map_step = m_size./map_size;
%%
populations_r = 1;
%populations_density=[20,10]; % people/km^2
%populations_positions=[20 20;70 70;]; % mapindex
means = 30;
populations_density = rand(1,means)*15;
populations_positions = unifrnd(1,50,means,2);
populations_risk = 10 * rand(map_size); % value
%%
road_density=[10,5,12,12]; % cars/km^@
road_positions=[25 25;40 25;10 20;40 30]; % mapindex
road_flag = [[ones(50,1)*25 [1:50]'];[ones(50,1)*24 [1:50]'];
             [[1:50]' ones(50,1)*11];[[1:50]' ones(50,1)*12];[[1:50]' ones(50,1)*13];
             [[1:50]' ones(50,1)*21];[[1:50]' ones(50,1)*22];
             [ones(50,1)*10 [1:50]' ];[ones(50,1)*9 [1:50]' ];[ones(50,1)*11 [1:50]' ];
             [ones(29,1)*37 [22:50]' ];[ones(29,1)*38 [22:50]' ];];
map_road_flag = zeros(size(map_size));
road_risk = zeros(map_size);
for i =1:size(road_flag,1)
    road_risk(road_flag(i,1),road_flag(i,2))=6;
    map_road_flag(road_flag(i,1),road_flag(i,2))=1;
end
%% update risk
for i=1:map_size(1)
    for j=1:map_size(2)
        pnum=length(populations_density);
        rnum=length(road_density);
        for k1=1:pnum
            r = sqrt(sum((([i,j]-populations_positions(k1,:)).*map_step).^2))/1000;
            if r<=populations_r
               populations_risk(i,j)=populations_risk(i,j)+exp(1-r^2)* populations_density(k1);
            end
        end
        for k2=1:rnum 
            r = sqrt(sum(([i,j]-road_positions(k2,:)).*map_step).^2)/1000;
            if r<=populations_r && map_road_flag(i,j)==1
                road_risk(i,j)=road_risk(i,j)+exp(1-r^2)* road_density(k2);
            end
        end
    end
end
%% buliding
miu = 3.0467;
sigma = 0.76023;
x = 1:1:140;
y = 1./(x.*sigma.*sqrt(2.*pi)).*exp(-(log(x)-miu).^2./(2*sigma^2));
figure();
xlabel('height[m]')
ylabel('P')


cumsum_y = [0 cumsum(y)+(1-sum(y))];
bulid_num = 200;
bulid_height = zeros(1,bulid_num);
for i =1:bulid_num
    index = find(rand()<cumsum_y,1,'first');
    bulid_height(i) = x(index-1); 
end
hold on
[n,xout]=hist(bulid_height,140);
bs=bar(xout,n/bulid_num);
for i =1:numel(bs)
    set(bs(i),'FaceColor','#9FB7A0');
end
plot(x,y,'black--','LineWidth',2);
bulid_position = round(unifrnd(1,50,bulid_num,2));
bulid_xyz = zeros(map_size);
for i =1:size(road_flag,1)
    bulid_xyz(road_flag(i,1),road_flag(i,2))=5;
end
for i =1:bulid_num
    if map_road_flag(bulid_position(i,1),bulid_position(i,2)) == 1
        continue;
    end
    bulid_xyz(bulid_position(i,1),bulid_position(i,2))=bulid_height(i);
    if rand<0.4
        if rand<0.5
            bulid_xyz(bulid_position(i,1)+1,bulid_position(i,2))=bulid_height(i);
        else
            bulid_xyz(bulid_position(i,1),bulid_position(i,2)+1)=bulid_height(i);
        end
    else
        if rand>0.7
            bulid_xyz(bulid_position(i,1)+1,bulid_position(i,2))=bulid_height(i);
            bulid_xyz(bulid_position(i,1),bulid_position(i,2)+1)=bulid_height(i);
            bulid_xyz(bulid_position(i,1)+1,bulid_position(i,2)+1)=bulid_height(i);
        else
        end
    end
    
end

figure()
hold on 
colormap('colorcube')
remove_empty_bars(bar3((bulid_xyz'),1));
%plot3(t(:,1),t(:,2),t(:,3),'Color','#A2AB28','LineStyle','-','Marker','o','LineWidth',1.5,'MarkerSize',5);
xlim([0 50])
ylim([0 50])
zlim([0 max(bulid_height)+10])
xlabel('X(100m/unit))');
ylabel('Y(100m/unit)');
zlabel('Z(m/unit)');
title('Buliding')

IOT_pos = [25,30,0;34,20,0;40,35 0;];
%%
figure()
subplot(1,2,1)
colormap('jet')
contourf(road_risk)
colorbar;
title('Road density')
subplot(1,2,2)
colormap('jet')
contourf(populations_risk)
title('Populations density')
colorbar;

