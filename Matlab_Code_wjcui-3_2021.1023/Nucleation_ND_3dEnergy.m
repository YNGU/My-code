%% ����matlab��ȡֵ���ٽ�Z������ճ����origin�л�ͼ
clc
clear
%% ����x��y��ֵ������Һ�ε��κ˹�ʽ
image = figure;  %�½�һ��ͼ��
%set(gcf,'unit','centimeters','position',[0 0 6 4]); %���û�����С����λ����
%set(gcf,'color',[1 1 1]) %���ñ�����ɫ��x/255=������ΧΪ0-1
r=0:0.05:10;  
h=0:0.05:10;  
V=1; %VΪ��Gv���趨ֵ
As=0.000000002; %AsΪ��s���趨ֵ
Ai=0.32; %AiΪ��i���趨ֵ
[X,Y]=meshgrid(r,h); 
for i=1:201;    %forѭ��ȥ������y>2x��ֵ��ǰ��x��y��Χ���˵Ļ�������ҲҪ�仯
    for j=1:201;
        if Y(i,j)>2.*X(i,j);  %(i,j)��ʾ�����еĶ�Ӧ��i�У���j�е�Ԫ��
           Y(i,j)=NaN;  %NaN��ʾȱʧֵ��û����ֵ
        end
    end
end
for i=1:401;    
    for j=1:401;
        if Y(i,j)<X(i,j);  %h<r�ķ�Χȥ��
           Y(i,j)=NaN;  
        end
    end
end
Gv=1.047.*V.*(Y.^2).*[(3*X)-Y];  %���ܹ��� �Թ�origin��ͼ��û����
GAs=(6.28*As*X).*Y;   %�����ܹ���
GAi=3.14*Ai.*(2.*Y.*X-Y.^2)   %(X.^2).*[1-((Y-X).^2)./(X.^2)];
%.^2���Ǿ�����ÿ��Ԫ�ص�ƽ��,�����ӦԪ�������.*,./Ҳһ��
Z=-Gv+GAs+GAi; %ͳһ��λ��nm��m֮��Ļ���
surf(X,Y,Z);
shading interp  %��ֵʹ��ͼ����ƽ��
title("ND ��G energy")
xlabel('Radius r'),ylabel('Height h'),zlabel('��G')
colormap(jet)
colorbar

%zticks([])��%��z��̶Ȳ���ʾ
hold on
savefig(image,'D:\++++++Data\Matlab Data\Exp_Data\ND_3dEnergy')


