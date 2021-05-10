function [ next ] = limit( current,next,Kinematic )
%UNTITLED 考虑加速和减速限制
%%Kinematic 最高速度m/s],最高旋转速度[rad/s],加速度[m/ss],旋转加速度[rad/ss]]
%% 速度上的限制
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
%% 角速度上的限制
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

