%% fitting ����excel����-1��Fick's Laws��
clc;
clear;
%%
data_x=xlsread('D:\++++++Data\Matlab Data\Exp_Data\GBT\Bi-content-Curve_fit\33.xlsx','sheet1','A1:A25');  %xlsread��excel������,excel���ڶ�ȡ·��
data_y=xlsread('D:\++++++Data\Matlab Data\Exp_Data\GBT\Bi-content-Curve_fit\33.xlsx','sheet1','B1:B25');  %����RPCԶ�̹��̵���ʧ�ܡ�ȥexcel������ļ�-ѡ��-����-com������-��add in�Ķ�ѡ��Ĺ���ȥ������
x=(data_x*0.38) %����x��nm����
y=data_y/255
%% �ֶ�����excel���ݣ�ѡ��ȷ�����ɵ�������-2
data_c=xlsread('D:\++++++Data\Matlab Data\Exp_Data\GBT\Bi-content-Curve_fit\33.xlsx',-1) 
%% ����txt��ʽ�ļ�
load x.txt
load y.txt

%% ��Ϸ���1
cftool %��������ϴ���

%% ��Ϸ���2 ʹ��interp1
plot(data_x,data_c,'o'); 
hold on;
x=min(data_x):.1:max(data_x);
c=interp1(data_x,data_c,x,'pchip'); % interp1�е�1��һ����
%������ѡ����������ѡ��ͬ����϶���ʽ����ʾ��ʽ�����ù�ʽ��ʾ����
plot(x,c);
hold off
%% generate code
%% ͨ����ʽ�õ��������,�����ֵ����ȥ
fit_x=[0:0.01:9.12];
fit_y=(-0.01018+0.9917)/2-((0.9917+0.01018)/2)*erf(0.6711*(fit_x-3.09));
plot(fit_x,fit_y);
fit_x1=fit_x';
fit_y1=fit_y';

