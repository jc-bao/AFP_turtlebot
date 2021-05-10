% function [ output_args ] = formation_avoidance( input_args )
%%%%一阶算法仿真 非完整型约束
%%% 2018-7-24
%%% 最初自己综合设计的方法，效果不好

    clear all;
    % 初始化 位置pose、速度V、加速度控制量control
       init_f=[-3 -6 0; %%%[x y th]
                -5 6 0;
                2 4 pi/4; 
                8 -3 -pi/4;
                2 0 pi/2];           
%     init_f=[0 0 0;
%             1 0 0;
%             2 0 0;
%             3 0 0;
%             4 0 0];
    pose_x=init_f(:,1);
    pose_y=init_f(:,2);
    pose_th=init_f(:,3);
    %% follower相对leader的位置
    delta_x=[-2 -6 -2 -6 0];   % 相对间隔误差   
    delta_y=[4 4 -4 -4 0];  %领航者与自己无误差

    fol_num=4;        
    N=5;             % 4follower and 1 leader
    countmax=1000;
    dt=0.1;
    
    K0=0;
    K_w=1;
    K_v=1;
    %%% 直线K1 K2都设在0.2左右 k3=0
    %%% 圆不考虑K2=0 K1，K3设为0.2左右
    K1=0.5;%%%位置偏差线速度调节
    K2=0;%%位置偏差角速度调节
    K3=0.3;%%朝向偏差角速度调节


    %% 通信拓扑图:1-4行为follower 最后一行为leader
    A=[0 1 1 1 1;     % a(ij)
       0 0 0 0 1;
       0 0 0 1 1;
       0 0 1 0 1;
       0 0 0 0 0];
% % %    %% 通信拓扑图:1-4行为follower 最后一行为leader
    A=[0 0 0 0 1;     % a(ij)
       0 0 0 0 1;
       0 0 0 0 1;
       0 0 0 0 1;
       0 0 0 0 0];

%        A=[0 1 1 1 3;     % a(ij)
%        0 0 0 0 3;
%        0 0 0 1 3;
%        0 0 1 0 3;
%        0 0 0 0 0];
    linear_v(:,1)=[0;0;0;0;0.5];
    angular_w(:,1)=[0;0;0;0;0.5];
    k=0;
    % 最高速度m/s],最高旋转速度[rad/s],加速度[m/ss],旋转加速度[rad/ss]]
    Kinematic=[1.5,toRadian(50.0),0.5,toRadian(50.0)];%% 运动限制
    flag=0;
    %% 开始循环 走顺时针圆周
%     figure;
    for count=1:countmax
        k=k+1;
        linear_v(N,k+1)=linear_v(N,k);%领航者速度不变
        angular_w(N,k+1)=angular_w(N,k);
        
%         if count==500
%             linear_v(N,k+1)=0.5;%领航者速度不变
%              angular_w(N,k+1)=0.5;
%         end
%         if count==1000
%             linear_v(N,k+1)=0.5;%领航者速度不变
%              angular_w(N,k+1)=0;
%         end
%         if (angular_w(N,k+1)~=0) %%%有转速
%             K2=0;
%             K3=0.2;
%         else
%             K2=0.2;
%             K3=0;
%         end
        
        for i=1:fol_num
            sum_delta_x=0;
            sum_delta_y=0;
            sum_delta_th=0;
            sum_delta_v=0;
            sum_delta_w=0;
            sum_ang=0;
            for j=1:N %%考虑邻居对它的影响
                sum_delta_x=sum_delta_x+A(i,j)*((pose_x(j,k)-pose_x(i,k))-(delta_x(j)-delta_x(i)));
                sum_delta_y=sum_delta_y+A(i,j)*((pose_y(j,k)-pose_y(i,k))-(delta_y(j)-delta_y(i)));
                sum_delta_th=sum_delta_th+A(i,j)*(pose_th(j,k)-pose_th(i,k));
                sum_delta_v=sum_delta_v+A(i,j)*linear_v(j,k);
                sum_delta_w=sum_delta_w+A(i,j)*angular_w(j,k);
            end
            sum_delta_v=sum_delta_v/sum(A(i,:));
            sum_delta_w=sum_delta_w/sum(A(i,:));
            delta_ang=atan2(sum_delta_y,sum_delta_x);
            delta_ang=seek_ang(delta_ang,pose_th(i,k));
%             if abs(delta_ang)<0.1
%                 delta_ang=0;
%             end

%%%正常版本
            angular_w(i,k+1)=K2*delta_ang+K0*angular_w(N,k)+K3*sum_delta_th+K_w*sum_delta_w;
            linear_v(i,k+1)=K1*(sum_delta_x*cos(pose_th(i,k))+sum_delta_y*sin(pose_th(i,k)))+K0*linear_v(N,k)+K_v*sum_delta_v;
            u_old=[linear_v(i,k);angular_w(i,k)];
            u=[linear_v(i,k+1);angular_w(i,k+1)];
            
            %%%加入速度限制
            u=limit(u_old,u,Kinematic);
          
            old_position=[pose_x(i,k);pose_y(i,k);pose_th(i,k)];
            new_position=motion(old_position,u,dt);
            pose_x(i,k+1)=new_position(1);
            pose_y(i,k+1)=new_position(2);
            pose_th(i,k+1)=new_position(3);
        end
        %% 更新领航者
        old_position=[pose_x(N,k);pose_y(N,k);pose_th(N,k)];
        u=[linear_v(N,k+1);angular_w(N,k+1)];
        new_position=motion(old_position,u,dt);
        pose_x(N,k+1)=new_position(1);
        pose_y(N,k+1)=new_position(2);
        pose_th(N,k+1)=new_position(3);

        tt_x(1:4,k)=pose_x(5,k);
        error_x(:,k)=tt_x(1:4,k)-pose_x(1:4,k)+(delta_x(1:4))';
        tt_y(1:4,k)=pose_y(5,k);
        error_y(:,k)=tt_y(1:4,k)-pose_y(1:4,k)+(delta_y(1:4))';
        
       %% ====Animation====
        hold off;
        ArrowLength=0.5;% 
        for j=1:N
            quiver(pose_x(j,k+1),pose_y(j,k+1),ArrowLength*cos(pose_th(j,k+1)),ArrowLength*sin(pose_th(j,k+1)),'*k');hold on;
            draw_circle (pose_x(j,k+1),pose_y(j,k+1),0.5);hold on;
        end
%         area=[-10 10 -10 10];
        area = compute_area(pose_x(N,k+1),pose_y(N,k+1),12);
        axis(area);
        grid on;
        drawnow;
    end
    %% 画误差图
    color='mgbkr'; %%%定义颜色标记
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
    hold on;
    grid on;
    xlabel('x');
    ylabel('y');
    legend('跟随者1','跟随者2','跟随者3','跟随者4','领航者');
    xlabel('x(m)');
    ylabel('y(m)');
    title('基于一致性的二阶编队算法');
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
% end

