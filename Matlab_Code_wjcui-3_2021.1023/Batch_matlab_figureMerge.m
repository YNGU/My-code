clc
clear

%% 批量的机械合并fig图像
for i=1:36    %fig的数量
    hf(i)=open(['D:\++++++Data\Matlab Data\Exp_Data\Bi liquid\simulatedTEM_moreAngle\02-21_241\', num2str(i),'.fig']);   %批量打开fig，输入路径
    fig(i)=get(hf(i), 'CurrentAxes');    %获取绘制的图像
end
image = figure;  %新建一个图像
set(gcf,'unit','centimeters','position',[0 0 18.3 20]); %设置画布大小，单位厘米
set(gcf,'color',[0.51 0.5 0.51]) %设置背景颜色，x/255=？，范围为0-0.1
for i = 1:36 
    subplot(6,6,i);   %包含第i个子图，第一个数为行数，第二个数字是列数
    axChildren = get(fig(i),'Children'); %获取绘制的图像
    copyobj(axChildren, gca);    %复制到当前的图窗中
    close(hf(i));   %关掉打开的单个图像
    colormap(gray);   %灰度显示图像
    axis image;    %固定图像原来的比列不变
    axis off;  %不显示坐标轴和示数
    hold on;
end
savefig(image,'D:\++++++Data\Matlab Data\Exp_Data\Bi liquid\simulatedTEM_moreAngle\02-21_241\img') %保存fig数据，tif要手动保存

