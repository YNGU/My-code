%% 先建立一个空矩阵-然后手动加入数据
clc
clear
marix_empty=zeros(70,1)
%% 找到统计图的拟合峰值对应的x轴
data=marix_empty
MAX=max(data)
MIN=min(data(find(data>0))) %找到大于零的最小值
logical=(data~=0)  %比较各个元素是否为0,然后进行分析：为0就显示为0，不是就显示为1
N=sum(logical(:)) %统计有多少个元素（数据）
%hist(data(find(data>0)),15) %hist(Y,X); X是将横坐标划分为几列Bar,统计Y在X这个区间划分下的个数
%hold on   %在以上基础上，再叠加下面的图
histfit(data(find(data>0)),15) %hisfit对直方图进行拟合 
obj=get(gca,'children') %导出对应你的图里面的线或者面的句柄Line，Bar
x=get(obj(1), 'xdata');  %提取拟合线obj(1)的x轴坐标,提取bar就是obj（2）
y=get(obj(1), 'ydata');  
[peak,index]=max(y)      %找到拟合线的最大值peak，和最大值对应的数组序号inedx 
Number_peak=x(:,index)    %再找到最大值对应的横坐标
Error=MAX-MIN  %最大最小值差距
%Zdisplay=[Number_peak, MAX, MIN, Error] %列出主要数据
%Zdisplay=Zdisplay'
ans=Number_peak