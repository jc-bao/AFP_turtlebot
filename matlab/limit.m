function [ next ] = limit( current,next,Kinematic )
%UNTITLED ���Ǽ��ٺͼ�������
%%Kinematic ����ٶ�m/s],�����ת�ٶ�[rad/s],���ٶ�[m/ss],��ת���ٶ�[rad/ss]]
%% �ٶ��ϵ�����
delta_v=next(1)-current(1);
if delta_v>=0
    next(1)=min(current(1)+delta_v,current(1)+Kinematic(3));
else
    next(1)=max(current(1)+delta_v,current(1)-Kinematic(3));
end
if next(1)>=0
    next(1)=min(next(1),Kinematic(1));
else
    next(1)=max(next(1),-Kinematic(1));
end
%% ���ٶ��ϵ�����
delta_w=next(2)-current(2);
if delta_w>=0
    next(2)=min(current(2)+delta_w,current(2)+Kinematic(4));
else
    next(2)=max(current(2)+delta_w,current(2)-Kinematic(4));
end
if next(2)>=0
    next(2)=min(next(2),Kinematic(2));
else
    next(2)=max(next(2),-Kinematic(2));
end


end

