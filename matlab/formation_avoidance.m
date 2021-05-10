% function [ output_args ] = formation_avoidance( input_args )
%%%%һ���㷨���� ��������Լ��
%%% 2018-7-24
%%% ����Լ��ۺ���Ƶķ�����Ч������

    clear all;
    % ��ʼ�� λ��pose���ٶ�V�����ٶȿ�����control
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
    %% follower���leader��λ��
    delta_x=[-2 -6 -2 -6 0];   % ��Լ�����   
    delta_y=[4 4 -4 -4 0];  %�캽�����Լ������

    fol_num=4;        
    N=5;             % 4follower and 1 leader
    countmax=1000;
    dt=0.1;
    
    K0=0;
    K_w=1;
    K_v=1;
    %%% ֱ��K1 K2������0.2���� k3=0
    %%% Բ������K2=0 K1��K3��Ϊ0.2����
    K1=0.5;%%%λ��ƫ�����ٶȵ���
    K2=0;%%λ��ƫ����ٶȵ���
    K3=0.3;%%����ƫ����ٶȵ���


    %% ͨ������ͼ:1-4��Ϊfollower ���һ��Ϊleader
    A=[0 1 1 1 1;     % a(ij)
       0 0 0 0 1;
       0 0 0 1 1;
       0 0 1 0 1;
       0 0 0 0 0];
% % %    %% ͨ������ͼ:1-4��Ϊfollower ���һ��Ϊleader
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
    % ����ٶ�m/s],�����ת�ٶ�[rad/s],���ٶ�[m/ss],��ת���ٶ�[rad/ss]]
    Kinematic=[1.5,toRadian(50.0),0.5,toRadian(50.0)];%% �˶�����
    flag=0;
    %% ��ʼѭ�� ��˳ʱ��Բ��
%     figure;
    for count=1:countmax
        k=k+1;
        linear_v(N,k+1)=linear_v(N,k);%�캽���ٶȲ���
        angular_w(N,k+1)=angular_w(N,k);
        
%         if count==500
%             linear_v(N,k+1)=0.5;%�캽���ٶȲ���
%              angular_w(N,k+1)=0.5;
%         end
%         if count==1000
%             linear_v(N,k+1)=0.5;%�캽���ٶȲ���
%              angular_w(N,k+1)=0;
%         end
%         if (angular_w(N,k+1)~=0) %%%��ת��
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
            for j=1:N %%�����ھӶ�����Ӱ��
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

%%%�����汾
            angular_w(i,k+1)=K2*delta_ang+K0*angular_w(N,k)+K3*sum_delta_th+K_w*sum_delta_w;
            linear_v(i,k+1)=K1*(sum_delta_x*cos(pose_th(i,k))+sum_delta_y*sin(pose_th(i,k)))+K0*linear_v(N,k)+K_v*sum_delta_v;
            u_old=[linear_v(i,k);angular_w(i,k)];
            u=[linear_v(i,k+1);angular_w(i,k+1)];
            
            %%%�����ٶ�����
            u=limit(u_old,u,Kinematic);
          
            old_position=[pose_x(i,k);pose_y(i,k);pose_th(i,k)];
            new_position=motion(old_position,u,dt);
            pose_x(i,k+1)=new_position(1);
            pose_y(i,k+1)=new_position(2);
            pose_th(i,k+1)=new_position(3);
        end
        %% �����캽��
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
    %% �����ͼ
    color='mgbkr'; %%%������ɫ���
    figure                                %   ������άƽ��ͼ  ����
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
    legend('������1','������2','������3','������4','�캽��');
    xlabel('x(m)');
    ylabel('y(m)');
    title('����һ���ԵĶ��ױ���㷨');
    %% �����ͼ
    cx=0:0.1:countmax/10;
    figure                                %   ������άƽ��ͼ  ����
    error=sqrt(error_x.^2+error_y.^2);
    for i=1:4
        plot(cx(1:countmax-1),error(i,1:countmax-1),color(1,i));
        hold on;
    end
    legend('������1','������2','������3','������4');
    xlabel('ʱ��(s)');
    ylabel('λ�����');
    title('�����������캽�ߵ��������');
% end

