%% ����matlab��ȡֵ���ٽ�Z������ճ����origin�л�ͼ
clc
clear
%% ����x��y��ֵ������Һ�ε��κ˹�ʽ
image = figure;  
d=0:0.00000000005:0.00000001;  %dΪ�����߿��
l=0:0.00000000005:0.00000001;  %lΪ�����߳���
V=1; 
As=0.000000002; 
Ai1=0.00000000032; %��������ֵӦ��ͬND�κ�ʱ�Ĳ�������һ��
Ai2=0.00000000032; %AiΪ��i2���趨ֵ
[X,Y]=meshgrid(d,l); 
%for i=1:201;    %forѭ��ȥ������y>2x��ֵ��ǰ��x��y��Χ���˵Ļ�������ҲҪ�仯
 %   for j=1:201;
  %      if Y(i,j)>2.*X(i,j);  %(i,j)��ʾ�����еĶ�Ӧ��i�У���j�е�Ԫ��
   %        Y(i,j)=NaN;  %NaN��ʾȱʧֵ��û����ֵ
    %    end
    %end
%end
Gv=0.785.*V.*(X.^2).*Y; 
GAs=(3.14*As*X).*Y;         
GAi=0.785*(Ai1+Ai2).*(X.^2);  %��������Ĺ���
Z=-Gv+GAs+GAi;
surf(X,Y,Z);
shading interp  %��ֵʹ��ͼ����ƽ��
title("ND ��G energy")
xlabel('Width d'),ylabel('Length l'),zlabel('��G/��Gv')
colormap(jet)
colorbar
hold on
savefig(image,'D:\++++++Data\Matlab Data\Exp_Data\NW_nucleation_jet')


