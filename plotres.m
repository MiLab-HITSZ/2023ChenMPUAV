function plotres(bulid_xyz,bulid_height,t,titlestr)
figure()
hold on 
colormap('colorcube')
remove_empty_bars(bar3(bulid_xyz));
plot3(t(:,1),t(:,2),t(:,3),'black-o','LineWidth',1.5);
xlim([0 100])
ylim([0 100])
zlim([0 max(bulid_height)+10])
xlabel('X(50m)');
ylabel('Y(50m)');
zlabel('Z(m)');
title(titlestr)
end