function [ feature_new ] = normalization( feature ,mi,ma,lower,upper)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
feature_new=lower+(upper-lower)*(feature-mi)/(ma-mi);
end

