function [ result ] = seek_ang( th1,th2 )
%UNTITLED4 ������
%   th1ΪҪת���Ľǡ�th2Ϊ��ǰλ�õĽǶ�
%  th1-th2��Ϊ��Բο�������ļн�

result=th1-th2;
if(th1>0)
    if(result>pi)
        result=result-2*pi;
    end
else
    if(result<-pi)
        result=result+2*pi;
    end
end


% if(th1*th2>=0)
%     result=th1-th2;
% else
%     if(th1>0)
%         result=th1-th2;
%         if(result>pi)
%             result=result-2*pi;
%         end
%     else
%         result=th1-th2;
%         if(result<-pi)
%             result=result+2*pi;
%         end
%     end
% end

end

