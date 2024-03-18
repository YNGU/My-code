%% clear
clear;
clc;
%% load mat file
load('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\Bi、Ge escape\30.mat');
%%

%% load position of txt files Fig.2
pos=load('D:\++++++Data\Matlab Data\Exp_Data\Test\16pos.txt'); %导入txt数据

folder='E:\'; %建一个文件夹；new_folder='D:\test'创建新文件夹


f=figure; %设置一个fig的窗口
f.PaperPosition=[1 1 7 14]; %Position由位置向量[left，bottom，width，height]组成，决定坐标轴位置
set(gca,'fontsize',9) %绘制字体大小
image=ReadDMFile('D:\++++++Data\Matlab Data\Exp_Data\Test\16pos.dm3'); %读DM文件，若只画点和线，除去这一句即可
image=double(image); 
imagesc(image); %将数组 image 中的数据显示为一个图像
axis image off %保持宽高比并取消坐标轴
ylim([0 2048]); %当前坐标轴或图表的 y 轴限制
xlim([0 2048]);
colormap(Dan_Mumford); %换色彩标尺
hold all %后续的绘图命令不仅保留之前的绘图结果，还记住当前使用的线型和颜色，从而新的绘图可以使用不同的颜色或线型以示区别。
%plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k'); %pos(:,1)为x,%MarkerEdgeColor：用于设置标记点的边框线条颜色;MarkerFaceColor：用于设置标记点的内部区域填充颜色;Markersize：用于设置标记点的大小
%plot(pos(:,2)+1,pos(:,1)-1,'o','markersize',3,'markerfacecolor','w','markeredgecolor','none');移动点对应的位置+1-1
x=pos(:,2); %读取二维矩阵的第一组数据，也就是x
y=pos(:,1); %读取二维矩阵的第二组数据，要是三维矩阵得括号里面要三个值，:表示
imagesc(image); %将数组 image 中的数据显示为一个图像
axis image off %保持宽高比并取消坐标轴
ylim([700 1530]); %当前坐标轴或图表的 y 轴限制
xlim([500 720]);
colormap(Dan_Mumford); %换色彩标尺
hold all

KNN=6; % K最邻近算法，是实现分类器中比较简单易懂的一种分类算法，找6个点

Idx=knnsearch(pos,pos,'K',KNN); %Idx = knnsearch(X,Y,Name,Value)；Idx = knnsearch(X,Y) 为Y中的每个查询点查找X中的最近邻居
Idx(:,1)=[]; %设置矩阵第一列为空值，0，也就是去掉近邻为自身的点 （:,1）idx的第一列数据

%Idx_angle=Idx-Idx; %这一步好像没啥用
for i=1:1:length(pos)  %for是循环条件，length是求向量的元素个数，也就是pos中有多少个点
    for j=1:1:KNN-1  %近邻点的数量，除了自己
        Idx_angle(i,j)=atan2(pos(Idx(i,j),1)-pos(i,1),pos(Idx(i,j),2)-pos(i,2)); %计算i和j之间的角度,atan2(x,y)第四象限
    end
end

range=0.5; %给的误差范围，如果键长太短或者太长需要改此参数
for i=1:1:length(pos) 
    % find the first vertex
    k1=find(abs(Idx_angle(i,:)-pi/2)<range); %找180度的点, abs返回一个数的绝对值
    if isempty(k1) %判断是否为空
        continue;
    end
    index1=Idx(i,k1(1));
    
    % find the second vertex
    k2=find(abs(Idx_angle(index1,:)-0.35)<range);  
    if isempty(k2)
        continue;
    end
    index2=Idx(index1,k2(1));
    
    % find the third vertex
    k3=find(abs(Idx_angle(index2,:)-(-pi/2))<range);
    if isempty(k3)
        continue;
    end
    index3=Idx(index2,k3(1));
    
    if mod(i,1)==0
        plot([pos(i,2) pos(index1,2) pos(index2,2) pos(index3,2) pos(i,2)],...
            [pos(i,1) pos(index1,1) pos(index2,1) pos(index3,1) pos(i,1)],...
            '-r','linewidth',0.7);  %plot画线
    end
    
end
%plot(x,y,'-c','linewidth',0.5) 
plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k');  %画点
print(f,'-dpng','-r300',['D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\1-5.png']); 
%%  --
%pos=load('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\Bi、Ge escape\29pos.txt');
x=pos(:,2); %读取二维矩阵的第一组数据，也就是x
y=pos(:,1); %读取二维矩阵的第二组数据，要是三维矩阵得括号里面要三个值，:表示

%% save to mat file
save('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\Bi、Ge escape\30.mat');

%% 以角度找第一个点
imagesc(image); %将数组 image 中的数据显示为一个图像
axis image off %保持宽高比并取消坐标轴
ylim([0 800]); %当前坐标轴或图表的 y 轴限制
xlim([400 600]);
colormap(inferno); %换色彩标尺
hold all

KNN=6;
Idx=knnsearch(pos,pos,'K',KNN);
Idx(:,1)=[];

Idx_angle=Idx-Idx;
for i=1:1:length(pos)
    for j=1:1:KNN-1
        Idx_angle(i,j)=atan2(pos(Idx(i,j),1)-pos(i,1),pos(Idx(i,j),2)-pos(i,2));
    end
end

range=0.2;
for i=1:1:length(pos)
    % find the first vertex
    k1=find(abs(Idx_angle(i,:)-pi/2)<range);
    if isempty(k1)
        continue;
    end
    index1=Idx(i,k1(1));
    
    % find the second vertex
    k2=find(abs(Idx_angle(index1,:)-0.35)<range);
    if isempty(k2)
        continue;
    end
    index2=Idx(index1,k2(1));
    
    % find the third vertex
    k3=find(abs(Idx_angle(index2,:)-(-pi/2))<range);
    if isempty(k3)
        continue;
    end
    index3=Idx(index2,k3(1));
    
    if mod(i,1)==0
        plot([pos(i,2) pos(index1,2) pos(index2,2) pos(index3,2) pos(i,2)],...
            [pos(i,1) pos(index1,1) pos(index2,1) pos(index3,1) pos(i,1)],...
            '-g','linewidth',0.5);
    end
    
end

%% 角度的分布直方图
index_ang=find(Idx_angle<150&Idx_angle>50) %找到小于170的数，有多少个
ang=Idx_angle(index_ang) 
%hist(Idx_angle(:),100); %找到角度分布直方图
hist(ang(:),100);

%% load position of txt files Fig.3
pos=load('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\001_surface\12-pos-2.txt'); %导入txt数据

f=figure; %设置一个fig的窗口
f.PaperPosition=[1 1 7 14]; %Position由位置向量[left，bottom，width，height]组成，决定坐标轴位置
set(gca,'fontsize',9) %绘制字体大小
%image=ReadDMFile('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\001_surface\1-inferno.dm3'); %读DM文件
image=double(image); 
imagesc(image); %将数组 image 中的数据显示为一个图像
axis image off %保持宽高比并取消坐标轴
ylim([0 1024]); %当前坐标轴或图表的 y 轴限制
xlim([300 500]);
colormap(inferno); %换色彩标尺
hold all %后续的绘图命令不仅保留之前的绘图结果，还记住当前使用的线型和颜色，从而新的绘图可以使用不同的颜色或线型以示区别。
%plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k'); %pos(:,1)为x,%MarkerEdgeColor：用于设置标记点的边框线条颜色;MarkerFaceColor：用于设置标记点的内部区域填充颜色;Markersize：用于设置标记点的大小
%plot(pos(:,2)+1,pos(:,1)-1,'o','markersize',3,'markerfacecolor','w','markeredgecolor','none');移动点对应的位置+1-1
x=pos(:,2); %读取二维矩阵的第一组数据，也就是x
y=pos(:,1); %读取二维矩阵的第二组数据，要是三维矩阵得括号里面要三个值，:表示
imagesc(image); %将数组 image 中的数据显示为一个图像
axis image off %保持宽高比并取消坐标轴
ylim([0 1024]); %当前坐标轴或图表的 y 轴限制
xlim([300 500]);
colormap(inferno); %换色彩标尺
hold all

KNN=6; % K最邻近算法，是实现分类器中比较简单易懂的一种分类算法
Idx=knnsearch(pos,pos,'K',KNN); %Idx = knnsearch(X,Y,Name,Value)；Idx = knnsearch(X,Y) 为Y中的每个查询点查找X中的最近邻居
Idx(:,1)=[]; %设置矩阵第一列为空值，0，也就是去掉近邻为自身的点

Idx_angle=Idx-Idx; %
for i=1:1:length(pos)  %for是循环条件，
    for j=1:1:KNN-1  %
        Idx_angle(i,j)=atan2(pos(Idx(i,j),1)-pos(i,1),pos(Idx(i,j),2)-pos(i,2)); %计算i和j之间的角度
    end
end

range=0.5; %给的误差范围，如果键长太短或者太长需要改此参数
for i=1:1:length(pos) 
    % find the first vertex
    k1=find(abs(Idx_angle(i,:)-pi/2)<range); %找180度的点
    if isempty(k1) %判断是否为空
        continue;
    end
    index1=Idx(i,k1(1));
    
    % find the second vertex
    k2=find(abs(Idx_angle(index1,:)-0.35)<range);
    if isempty(k2)
        continue;
    end
    index2=Idx(index1,k2(1));
    
    % find the third vertex
    k3=find(abs(Idx_angle(index2,:)-(-pi/2))<range);
    if isempty(k3)
        continue;
    end
    index3=Idx(index2,k3(1));
    
    if mod(i,1)==0
        plot([pos(i,2) pos(index1,2) pos(index2,2) pos(index3,2) pos(i,2)],...
            [pos(i,1) pos(index1,1) pos(index2,1) pos(index3,1) pos(i,1)],...
            '-k','linewidth',1);
    end
    
end
%plot(x,y,'-c','linewidth',0.5) 
plot(pos(:,2),pos(:,1),'o','markersize',3,'markerfacecolor','g','markeredgecolor','k');
print(f,'-dpng','-r300',['D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\001_surface\12plus-pos.png']);
%% pos

