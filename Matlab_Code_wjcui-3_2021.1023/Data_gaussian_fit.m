%% �Ƚ���һ���վ���-Ȼ���ֶ���������
clc
clear
marix_empty=zeros(70,1)
%% �ҵ�ͳ��ͼ����Ϸ�ֵ��Ӧ��x��
data=marix_empty
MAX=max(data)
MIN=min(data(find(data>0))) %�ҵ����������Сֵ
logical=(data~=0)  %�Ƚϸ���Ԫ���Ƿ�Ϊ0,Ȼ����з�����Ϊ0����ʾΪ0�����Ǿ���ʾΪ1
N=sum(logical(:)) %ͳ���ж��ٸ�Ԫ�أ����ݣ�
%hist(data(find(data>0)),15) %hist(Y,X); X�ǽ������껮��Ϊ����Bar,ͳ��Y��X������仮���µĸ���
%hold on   %�����ϻ����ϣ��ٵ��������ͼ
histfit(data(find(data>0)),15) %hisfit��ֱ��ͼ������� 
obj=get(gca,'children') %������Ӧ���ͼ������߻�����ľ��Line��Bar
x=get(obj(1), 'xdata');  %��ȡ�����obj(1)��x������,��ȡbar����obj��2��
y=get(obj(1), 'ydata');  
[peak,index]=max(y)      %�ҵ�����ߵ����ֵpeak�������ֵ��Ӧ���������inedx 
Number_peak=x(:,index)    %���ҵ����ֵ��Ӧ�ĺ�����
Error=MAX-MIN  %�����Сֵ���
%Zdisplay=[Number_peak, MAX, MIN, Error] %�г���Ҫ����
%Zdisplay=Zdisplay'
ans=Number_peak