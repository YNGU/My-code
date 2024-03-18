%% fitting 导入excel数据-1（Fick's Laws）
clc;
clear;
%%
data_x=xlsread('D:\++++++Data\Matlab Data\Exp_Data\GBT\Bi-content-Curve_fit\33.xlsx','sheet1','A1:A25');  %xlsread读excel的数据,excel放在读取路径
data_y=xlsread('D:\++++++Data\Matlab Data\Exp_Data\GBT\Bi-content-Curve_fit\33.xlsx','sheet1','B1:B25');  %报错“RPC远程过程调用失败”去excel里面的文件-选项-管理-com加载项-带add in的都选项的勾都去掉就行
x=(data_x*0.38) %计算x的nm距离
y=data_y/255
%% 手动导入excel数据，选择并确定即可导入数据-2
data_c=xlsread('D:\++++++Data\Matlab Data\Exp_Data\GBT\Bi-content-Curve_fit\33.xlsx',-1) 
%% 导入txt格式文件
load x.txt
load y.txt

%% 拟合方法1
cftool %打开曲线拟合窗口

%% 拟合方法2 使用interp1
plot(data_x,data_c,'o'); 
hold on;
x=min(data_x):.1:max(data_x);
c=interp1(data_x,data_c,x,'pchip'); % interp1中的1是一不是
%工具中选择基本拟合来选择不同的拟合多项式，显示公式可以让公式显示出来
plot(x,c);
hold off
%% generate code
%% 通过公式得到拟合曲线,把拟合值带进去
fit_x=[0:0.01:9.12];
fit_y=(-0.01018+0.9917)/2-((0.9917+0.01018)/2)*erf(0.6711*(fit_x-3.09));
plot(fit_x,fit_y);
fit_x1=fit_x';
fit_y1=fit_y';

