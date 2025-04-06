close all
clear
clc
%%
load('traj.mat')
load('PredSteps.mat')

dimen   = size(traj);
max_idx = dimen(2);
Nv      = dimen(3);

R = 0.3;
r = R*0.2;
Lw = R*0.6;
hw = R*0.3;

str = {'r', 'y','c',' m','b', [0.6350 0.0780 0.1840], [0 0.4470 0.7410], [0.4940 0.1840 0.5560]};
% set(gcf, 'Position',  [2400, 400, 800, 800]);
for snap = 1:2:max_idx
    clf(figure(snap));
    set(gca, 'position', [0.06 0.06 0.92 0.92], 'FontSize', 10);
    t = tiledlayout(1,1,'Padding','tight');
    t.Units = 'inches';
    t.OuterPosition = [0.25 0.25 3 3];
    nexttile;

    for i = 1:Nv
        circle(traj(1,snap,i), traj(2,snap,i), R, str{i});
        hold on
        circle(traj(1,snap,i)+(R-r)*cos(traj(3,snap,i)),traj(2,snap,i)+(R-r)*sin(traj(3,snap,i)),r,'w');
        hold on
        plot(PredSteps(1,:,snap,i), PredSteps(2,:,snap,i),'LineWidth', 1, 'Color','g')
        hold on
        Vehicle([traj(1,snap,i)+(R-0.9*hw)*cos(traj(3,snap,i)+pi/2),traj(2,snap,i)+(R-0.9*hw)*sin(traj(3,snap,i)+pi/2)],Lw,hw,traj(3,snap,i),'black');
        hold on
        Vehicle([traj(1,snap,i)+(R-0.9*hw)*cos(traj(3,snap,i)-pi/2),traj(2,snap,i)+(R-0.9*hw)*sin(traj(3,snap,i)-pi/2)],Lw,hw,traj(3,snap,i),'black');    
        hold on
        plot(traj(1,end,i), traj(2,end,i), 'k+','LineWidth',2)
    end
    xlim([2 10]);
    ylim([2 10]);
    %exportgraphics(t,append('figs1/snap',num2str(snap-1),'.pdf'), 'ContentType','vector');
    pause(0.5)   
end

xplot = zeros(1,max_idx-1);
for k = 2:max_idx-1
xplot(k) = norm(PredSteps(:,1,k,1) - traj(1:2,k,1));
end

figure(100)
hold on
plot(1:max_idx-1, xplot);