%% 循环填充四面体
%第一个点的位置
clc
clear all
fig=openfig('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\test\vdw gap-octahedrons_POS.fig')
pos=load('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\test\vdw gap-octahedrons_POS.txt'); %导入txt数据
y_min=min(pos(:,2));    %数组索引必须为正整数或逻辑值出现：用clear all解决
range=5;    %设置一个误差范围 %ctrl+c强制执行：可以中止继续输入语句
A=find(pos(:,2)<y_min+range);  %find找到y极小值对应的行数
x_min=min(pos(A,1));    %找到y值最小对应的x值最小值
[row,column]=find(pos(:,1)==x_min);  %find找到值对应的行和列
y=pos(row,2);   %最小x值对应的y值
pos_first=[x_min,y] %第一个点的位置

%找第一个点的四个近邻点并画四边形
KNN=4; % K最邻近算法，找带自身外的四个点
Idx=knnsearch(pos,pos_first,'K',KNN); %Idx = knnsearch(X,Y,Name,Value) 为Y中的每个查询点查找X中的最近邻
four_spots=pos(Idx,:); %找出4个近邻点的坐标
[~,l]=sort(four_spots(:,1))   %对四个点之和进行排序
four_spots=four_spots(l,:)  
M=[four_spots(1),four_spots(2),four_spots(4),four_spots(3)]; %按照1243的顺序排列x坐标
N=[four_spots(5),four_spots(6),four_spots(8),four_spots(7)];  %按照1243对应的顺序y坐标
patch(M,N,'r')  %将四个点按顺序填充颜色
alpha(0.3)  %alpha 函数设置当前坐标区范围内所有图像、填充或曲面对象的透明度。
hold all
%axis image

%找第二个顶点
times=40;      %循环次数
for i=1:times;   %for循环对后面的语句进行循环9次的操作
pos(row,:)=[];  %在pos中去除上面已经画过的第一个点
y2_min=min(pos(:,2));    %数组索引必须为正整数或逻辑值出现：用clear all解决
range=5;    %设置一个误差范围 %ctrl+c强制执行：可以中止继续输入语句
A=find(pos(:,2)<y_min+range);  %find找到y极小值对应的行数
x2_min=min(pos(A,1));    %找到y值最小对应的x值最小值
[row,column]=find(pos(:,1)==x2_min);  %find找到值对应的行和列
y2=pos(row,2);   %最小x值对应的y值
pos_second=[x2_min,y2] %第二个点的位置

%找到第二个点的5个近邻点，并删除多余的点，做四边形填充
KNN=5; % K最邻近算法，找带自身外的五个点
Idx=knnsearch(pos,pos_second,'K',KNN); %Idx = knnsearch(X,Y,Name,Value) 为Y中的每个查询点查找X中的最近邻
five_spots=pos(Idx,:);    %找出4个近邻点的坐标
extra_spot=find(five_spots(:,1)<x2_min); %找出多余点的x坐标
five_spots(extra_spot,:)=[];   %去掉多余点 
four_spots=five_spots    %得到剩下的四个点
[~,l]=sort(four_spots(:,1))   %对四个点之和进行排序
four_spots=four_spots(l,:)  
M=[four_spots(1),four_spots(2),four_spots(4),four_spots(3)]; %按照1243的顺序排列x坐标
N=[four_spots(5),four_spots(6),four_spots(8),four_spots(7)];  %按照1243对应的顺序y坐标
patch(M,N,'r')  %将四个点按顺序填充颜色
hold all
alpha(0.3)  %透明度
end    %对循环次数语句end

%% 从不同方向找四个近邻点，并计算距离与角度作文填充值
clc
clear all
fig=openfig('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\test\vdw gap-octahedrons_POS.fig')
pos=load('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\test\vdw gap-octahedrons_POS.txt'); %导入txt数据
y_min=min(pos(:,2));    %数组索引必须为正整数或逻辑值出现：用clear all解决
range=5;    %设置一个误差范围 %ctrl+c强制执行：可以中止继续输入语句
A=find(pos(:,2)<y_min+range);  %find找到y极小值对应的行数
x_min=min(pos(A,1));    %找到y值最小对应的x值最小值
[row,column]=find(pos(:,1)==x_min);  %find找到值对应的行和列
y=pos(row,2);   %最小x值对应的y值
pos_first=[x_min,y] %第一个点的位置。

KNN=4; %找第一个点
Idx=knnsearch(pos,pos_first,'K',KNN);
four_spots=pos(Idx,:); %四个点坐标
for i=1:1:length(four_spots);
    A=four_spots(i,:)-pos_first  %四个点与第一个点的差值    
end    
B=[atan2(A(1),A(5));atan2(A(2),A(6));atan2(A(3),A(7));atan2(A(4),A(8))] %找到四个点的角度
range_angle=0.2
find(B,pi/2-2<B<pi/2+0.2)









