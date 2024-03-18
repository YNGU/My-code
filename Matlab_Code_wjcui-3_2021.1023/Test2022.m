%% ZnO-Al2O3_Xupengyu_2022.0214 矢量关系
h=2
k=0
l=0
h1=0.5*h-0.5*l;
k1=0.5*k-0.5*h;
l1=h+k+l;
ans=[h1,k1,l1]

%%
clc
clear
data_per=xlsread('D:\++++++Data\1\per.xlsx','sheet1','B1:B770');  %xlsread读excel的数据,excel放在读取路径
data_def=xlsread('D:\++++++Data\1\vsb.xlsx','sheet1','B1:B770');  
data=data_per-data_def
%% a版本稳定性不高，b版稳定性比较高---appdesigner打开app设计
a=1+1;   
b=1.*10;    %正常的乘除都加.符号,还有平方根,次方根都要加
c=b./a
%% 零基础matlab系列教程 P1-P2____number
%command window
 %help+function / 右上角search找想要的函数
cos(((1+2+3+4)^3/5)^0.5);
sin(sqrt(pi))+log(tan(1));
2^(3.5*1.7);
exp(sin(10));
cos(pi):
sin(ans):
  %Variables变量
  A=10; 
  a=20;
  A2=10;   %variable/data type, who看现有哪些变量,whos看现有的变量和types/
  %ans,i,j complex number,inf,eps 2.2204e-016,NaN not a
  %number,pi,都已经有含义的variables,叫iskeywords
  cos='This string.'; %这里的cos变成了variables,不是cos函数了
  cos(8);  %这里的cos(8)是指cos中第八个字母;variables比built-in function的优先级更高
 %format long / shortE/longE/bank/hex/rat
%% 零基础matlab系列教程 P1-P2____array
a=[1 2 3 4]; %行向量 row vector
b=[1; 2; 3; 4]; %列向量 column vector
c=a*b;
d=b*a;
A=[1, 21, 6; 5, 17, 9; 31, 2, 7];
a(3);
b(4);  %array indexing 找矩阵里面的数字
A(1, 2);   %找矩阵的位置数值
A(2);
A(8);
A([1 3 5]);
A([1 3; 1 3]);
A([1 3], [1 3]); %前面的方括号是找第一行与第三行，后面的第一列和第三列。
A(1,2)=16;
A([1,2],[2,3])=[0 0; 0 0]
A(3,:)={}; %让一行全部为0
F=[A B]; %增广矩阵
A=[1 2 3; 4 5 4; 9 8 7]
B=[3 3 3; 2 4 9; 1 3 1]   %+，-，*，/  矩阵的加减乘除算法,.*点乘是对应位置相乘,A/B=A*inv(B);A./B=对应的位置数值相除；'转置矩阵
a=2
A^a
A.^a
%some special Matrix
linspace(0 13 5)
eye(3)
zeros(4,2)
ones(4,2)
diag([2 3 4]) %对角线矩阵
rand(3)
max(A)
max(max(A))
min(A)
sum(A)
sum(sum(A))
mean(A)
sort(A)
sortrows(A)
size(A)
length(A)
find(A==6)
%% load position of txt files Fig.2
pos=load('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\Te1-Sb-Te1_octahedrons\Te1-Sb-Te1_octahedrons-Pos-28.txt'); %导入txt数据

folder='E:\'; %建一个文件夹；new_folder='D:\test'创建新文件夹


f=figure; %设置一个fig的窗口
f.PaperPosition=[1 1 7 14]; %Position由位置向量[left，bottom，width，height]组成，决定坐标轴位置
set(gca,'fontsize',9) %绘制字体大小
image=ReadDMFile('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\1607 2.9 Mx STEM Diffraction HAADF 20211107 Nano-G.dm3'); %读DM文件，若只画点和线，除去这一句即可
image=double(image); 
%imagesc(image); %将数组 image 中的数据显示为一个图像
%axis image off %保持宽高比并取消坐标轴
%ylim([0 2048]); %当前坐标轴或图表的 y 轴限制
%xlim([0 2048]);
%colormap(inferno); %换色彩标尺
%hold all %后续的绘图命令不仅保留之前的绘图结果，还记住当前使用的线型和颜色，从而新的绘图可以使用不同的颜色或线型以示区别。
%plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k'); %pos(:,1)为x,%MarkerEdgeColor：用于设置标记点的边框线条颜色;MarkerFaceColor：用于设置标记点的内部区域填充颜色;Markersize：用于设置标记点的大小
%plot(pos(:,2)+1,pos(:,1)-1,'o','markersize',3,'markerfacecolor','w','markeredgecolor','none');移动点对应的位置+1-1
x=pos(:,2); %读取二维矩阵的第一组数据，也就是x
y=pos(:,1); %读取二维矩阵的第二组数据，要是三维矩阵得括号里面要三个值，:表示
imagesc(image); %将数组 image 中的数据显示为一个图像
axis image off %保持宽高比并取消坐标轴
ylim([0 2048]); %当前坐标轴或图表的 y 轴限制
xlim([0 2048]);
colormap(inferno); %换色彩标尺
hold all

KNN=6; % K最邻近算法，是实现分类器中比较简单易懂的一种分类算法，找6个点

Idx=knnsearch(pos,pos,'K',KNN); %Idx = knnsearch(X,Y,Name,Value)；Idx = knnsearch(X,Y) 为Y中的每个查询点查找X中的最近邻居
Idx(:,1)=[]; %设置矩阵第一列为空值，0，也就是去掉近邻为自身的点 （:,1）idx的第一列数据

%Idx_angle=Idx-Idx; %这一步好像没啥用
for i=1:1:length(pos)  %for是循环条件，length是求向量的元素个数，设定i是pos中的第几个位点
    for j=1:1:KNN-1  %近邻点的数量，除了自己，设定j为KNN中的点
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
            '-r','linewidth',0.7);  %plot画线，给找出的八面体四个点连线
    end
    if mod(i,1)==0
        c=rand(28); %添加颜色变量
        patch([pos(i,2) pos(index1,2) pos(index2,2) pos(index3,2) pos(i,2)],...
            [pos(i,1) pos(index1,1) pos(index2,1) pos(index3,1) pos(i,1)],...
            c);  %patch填充四个点的平行四边形为红色
        colormap(inferno)
    end
    
end
plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k');  %给原子位置画点
print(f,'-dpng','-r300',['D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\1-5-test.png']); 
%% 填充四个点的四面体
patch(pos(:,1),pos(:,2),'red');
imagesc(imgae)
hold all;


%% for循环语句
clear;
sum=0;
for i=1:100;
    sum=sum+i;
end
sum=sum %简单的for循环语句

clear;
sum=0;
for i=1:2:100;   %2是设置的步长
    sum=sum+i;
end
sum=sum    %加上步长的for循环语句

clear;
A=rand(1,4);
for i=A
    i
end            %对向量的遍历：指沿着某条搜索路线，依次对树中每个结点均做一次且仅做一次访问

clear;
for i=1:1:5
    for j=1:1:3
        A=i*j
    end
end

%% if 循环语句
clear
a=4
b=2
c=10
if(a<b)
    a=b
    if(a<c)
    a=c
    end
end 


clear;
a=12;
b=15;
c=18;
if(a<b)
    b=a
    if(b>c)
        b=c
    else
        b=c-a;
    end
end

%% fill
clear all
x=[0 1 1 0 0];
y=[0 0 1 1 0];
fill(x,y,'r')
axis([0 3 0 3])
%% 
x=[0 1 1 0]
y=[0 0 1 1]
patch(x,y,'red')

%% 
dis=load('C:\Users\Lenovo\Desktop\dis.txt');

 
  
  