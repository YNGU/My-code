clc
clear

%% �����Ļ�е�ϲ�figͼ��
for i=1:36    %fig������
    hf(i)=open(['D:\++++++Data\Matlab Data\Exp_Data\Bi liquid\simulatedTEM_moreAngle\02-21_241\', num2str(i),'.fig']);   %������fig������·��
    fig(i)=get(hf(i), 'CurrentAxes');    %��ȡ���Ƶ�ͼ��
end
image = figure;  %�½�һ��ͼ��
set(gcf,'unit','centimeters','position',[0 0 18.3 20]); %���û�����С����λ����
set(gcf,'color',[0.51 0.5 0.51]) %���ñ�����ɫ��x/255=������ΧΪ0-0.1
for i = 1:36 
    subplot(6,6,i);   %������i����ͼ����һ����Ϊ�������ڶ�������������
    axChildren = get(fig(i),'Children'); %��ȡ���Ƶ�ͼ��
    copyobj(axChildren, gca);    %���Ƶ���ǰ��ͼ����
    close(hf(i));   %�ص��򿪵ĵ���ͼ��
    colormap(gray);   %�Ҷ���ʾͼ��
    axis image;    %�̶�ͼ��ԭ���ı��в���
    axis off;  %����ʾ�������ʾ��
    hold on;
end
savefig(image,'D:\++++++Data\Matlab Data\Exp_Data\Bi liquid\simulatedTEM_moreAngle\02-21_241\img') %����fig���ݣ�tifҪ�ֶ�����

