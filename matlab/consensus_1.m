function [ output_args ] = consensus_1( input_args )
%%%验证一阶一致性算法 完整型约束 有避障
%%% 2018-7-25
close all;
addpath(genpath('./subfun'))
fol_num=4;        
N=5;             % 4follower and 1 leader
countmax=200;
dt=0.1;
gama=2;
beta=100;
K0=1;
% x最高速度m/s],y最高旋转速度[rad/s],x最高加速度[m/ss],y最高加速度[rad/ss]]
Kinematic=[1;1;0.2;0.2];
%% 1-4行为follower 最后一行为leader
A=[0 1 1 1 1;     % a(ij)
   0 0 0 0 1;
   0 0 0 1 1;
   0 0 1 0 1;
   0 0 0 0 0];
 %% 初始化 位置pose、速度V、加速度控制量control
    init_f=[-3 -3 0; %%%[x y th]
                -4 2 0;
                2 4 pi/4; 
                8 -3 -pi/4;
                0 0 pi/2];           
    pose_x=init_f(:,1);
    pose_y=init_f(:,2);
    pose_th=init_f(:,3);
    %% follower相对leader的位置
    delta_x=[-2 -6 -2 -6 0];   % 相对间隔误差   
    delta_y=[4 4 -4 -4 0];  %领航者与自己无误差
    V_x(:,1)=[0;0;0;0;1];
    V_y(:,1)=[0;0;0;0;1]; %%%leader在y方向的初始速度为1m/s
    k=0;
    %% 开始循环 走顺时针圆周
    for count=1:countmax
        k=k+1;
%         %%%做直线
%         V_x(N,k+1)=V_x(N,k);
%         V_y(N,k+1)=V_y(N,k);
        %%%做圆周
        V_x(N,k+1)=cos(k*dt);
        V_y(N,k+1)=sin(k*dt);
        for i=1:fol_num        
            sum_delta_x=0;
            sum_delta_y=0;
            for j=1:N %%考虑邻居对它的影响
                sum_delta_x=sum_delta_x+A(i,j)*((pose_x(j,k)-pose_x(i,k))-(delta_x(j)-delta_x(i)));
                sum_delta_y=sum_delta_y+A(i,j)*((pose_y(j,k)-pose_y(i,k))-(delta_y(j)-delta_y(i)));   
            end
           %%%考虑冲突避免加上斥力
            kk=0;
            for j=1:N
                if j~=i
                    kk=kk+1;
                    obs_pose(kk,1)=pose_x(j,k);
                    obs_pose(kk,2)=pose_y(j,k);
                end
            end
            repulsion=compute_repulsion([pose_x(i,k),pose_y(i,k)],obs_pose,2);        
            %%%%%
            V_x(i,k+1)=K0*V_x(N,k)+gama*sum_delta_x+beta*repulsion(1);
            V_y(i,k+1)=K0*V_y(N,k)+gama*sum_delta_y+beta*repulsion(2);
            out=confine([V_x(i,k) V_y(i,k)],[V_x(i,k+1) V_y(i,k+1)],Kinematic);
            V_x(i,k+1)=out(1);
            V_y(i,k+1)=out(2);
        end
        for i=1:N
            pose_x(i,k+1)=pose_x(i,k)+dt*V_x(i,k+1);
            pose_y(i,k+1)=pose_y(i,k)+dt*V_y(i,k+1);
            pose_th(j,k+1)=0;
        end
        tt_x(1:4,k)=pose_x(5,k);
        error_x(:,k)=tt_x(1:4,k)-pose_x(1:4,k)+(delta_x(1:4))';
        tt_y(1:4,k)=pose_y(5,k);
        error_y(:,k)=tt_y(1:4,k)-pose_y(1:4,k)+(delta_y(1:4))';
        %% ====Animation====
        hold off;
        ArrowLength=0.1;% 
        for j=1:N
            quiver(pose_x(j,k+1),pose_y(j,k+1),ArrowLength*cos(pose_th(j,k+1)),ArrowLength*sin(pose_th(j,k+1)),'*k');hold on;
            draw_circle (pose_x(j,k+1),pose_y(j,k+1),0.5);hold on;
        end
        area=[-10 10 -10 10];
        axis(area);
        grid on;
        drawnow;       
    end

    color='mgbkr'; %%%定义颜色标记
    %% 画图
    figure                                %   生成三维平面图  连续
    for i=1:N
        plot(pose_x(i,:),pose_y(i,:),color(1,i),'LineWidth',2);
        hold on
    end
    for i=1:N-1
        plot(pose_x(i,1),pose_y(i,1),'bp','color',color(1,i),'LineWidth',1);
        hold on
    end
    plot(pose_x(N,1),pose_y(N,1),'*','color',color(1,N),'LineWidth',1);
    hold on
    for i=1:N-1
        plot(pose_x(i,countmax),pose_y(i,countmax),'m^','color',color(1,i),'LineWidth',2);
        hold on
    end
    plot(pose_x(N,countmax),pose_y(N,countmax),'o','color',color(1,N),'LineWidth',2);
    hold on

    grid on;
    xlabel('x');
    ylabel('y');
    legend('跟随者1','跟随者2','跟随者3','跟随者4','领航者');
    xlabel('x(m)');
    ylabel('y(m)');
    title('基于一致性的一阶编队算法');
    %% 画误差图
    cx=0:0.1:countmax/10;
    figure                                %   生成三维平面图  连续
    error=sqrt(error_x.^2+error_y.^2);
    for i=1:4
        plot(cx(1:countmax-1),error(i,1:countmax-1),color(1,i));
        hold on;
    end
    legend('跟随者1','跟随者2','跟随者3','跟随者4');
    xlabel('时间(s)');
    ylabel('位置误差');
    title('各跟随者与领航者的误差曲线');
end

function [ next] = confine(current,next,Kinematic)
%%%current=[v_x v_y];
%%%%Kinematic=[ x最高速度m/s],y最高速度[m/s],x最高加速度[m/ss],y最高加速度[m/ss]]
%%%Kinematic=[1;1;0.5;0.5];
%% 速度x上的限制
delta_x=next(1)-current(1);
if delta_x>=0
    next(1)=min(current(1)+delta_x,current(1)+Kinematic(3));
else
    next(1)=max(current(1)+delta_x,current(1)-Kinematic(3));
end
if next(1)>=0
    next(1)=min(next(1),Kinematic(1));
else
    next(1)=max(next(1),-Kinematic(1));
end
%% 速度y上的限制
delta_y=next(2)-current(2);
if delta_y>=0
    next(2)=min(current(2)+delta_y,current(2)+Kinematic(4));
else
    next(2)=max(current(2)+delta_y,current(2)-Kinematic(4));
end
if next(2)>=0
    next(2)=min(next(2),Kinematic(2));
else
    next(2)=max(next(2),-Kinematic(2));
end
end

