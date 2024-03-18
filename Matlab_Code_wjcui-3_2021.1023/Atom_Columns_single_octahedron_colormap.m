%% ѭ�����������
%��һ�����λ��
clc
clear all
fig=openfig('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\test\vdw gap-octahedrons_POS.fig')
pos=load('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\test\vdw gap-octahedrons_POS.txt'); %����txt����
y_min=min(pos(:,2));    %������������Ϊ���������߼�ֵ���֣���clear all���
range=5;    %����һ����Χ %ctrl+cǿ��ִ�У�������ֹ�����������
A=find(pos(:,2)<y_min+range);  %find�ҵ�y��Сֵ��Ӧ������
x_min=min(pos(A,1));    %�ҵ�yֵ��С��Ӧ��xֵ��Сֵ
[row,column]=find(pos(:,1)==x_min);  %find�ҵ�ֵ��Ӧ���к���
y=pos(row,2);   %��Сxֵ��Ӧ��yֵ
pos_first=[x_min,y] %��һ�����λ��

%�ҵ�һ������ĸ����ڵ㲢���ı���
KNN=4; % K���ڽ��㷨���Ҵ���������ĸ���
Idx=knnsearch(pos,pos_first,'K',KNN); %Idx = knnsearch(X,Y,Name,Value) ΪY�е�ÿ����ѯ�����X�е������
four_spots=pos(Idx,:); %�ҳ�4�����ڵ������
[~,l]=sort(four_spots(:,1))   %���ĸ���֮�ͽ�������
four_spots=four_spots(l,:)  
M=[four_spots(1),four_spots(2),four_spots(4),four_spots(3)]; %����1243��˳������x����
N=[four_spots(5),four_spots(6),four_spots(8),four_spots(7)];  %����1243��Ӧ��˳��y����
patch(M,N,'r')  %���ĸ��㰴˳�������ɫ
alpha(0.3)  %alpha �������õ�ǰ��������Χ������ͼ��������������͸���ȡ�
hold all
%axis image

%�ҵڶ�������
times=40;      %ѭ������
for i=1:times;   %forѭ���Ժ����������ѭ��9�εĲ���
pos(row,:)=[];  %��pos��ȥ�������Ѿ������ĵ�һ����
y2_min=min(pos(:,2));    %������������Ϊ���������߼�ֵ���֣���clear all���
range=5;    %����һ����Χ %ctrl+cǿ��ִ�У�������ֹ�����������
A=find(pos(:,2)<y_min+range);  %find�ҵ�y��Сֵ��Ӧ������
x2_min=min(pos(A,1));    %�ҵ�yֵ��С��Ӧ��xֵ��Сֵ
[row,column]=find(pos(:,1)==x2_min);  %find�ҵ�ֵ��Ӧ���к���
y2=pos(row,2);   %��Сxֵ��Ӧ��yֵ
pos_second=[x2_min,y2] %�ڶ������λ��

%�ҵ��ڶ������5�����ڵ㣬��ɾ������ĵ㣬���ı������
KNN=5; % K���ڽ��㷨���Ҵ�������������
Idx=knnsearch(pos,pos_second,'K',KNN); %Idx = knnsearch(X,Y,Name,Value) ΪY�е�ÿ����ѯ�����X�е������
five_spots=pos(Idx,:);    %�ҳ�4�����ڵ������
extra_spot=find(five_spots(:,1)<x2_min); %�ҳ�������x����
five_spots(extra_spot,:)=[];   %ȥ������� 
four_spots=five_spots    %�õ�ʣ�µ��ĸ���
[~,l]=sort(four_spots(:,1))   %���ĸ���֮�ͽ�������
four_spots=four_spots(l,:)  
M=[four_spots(1),four_spots(2),four_spots(4),four_spots(3)]; %����1243��˳������x����
N=[four_spots(5),four_spots(6),four_spots(8),four_spots(7)];  %����1243��Ӧ��˳��y����
patch(M,N,'r')  %���ĸ��㰴˳�������ɫ
hold all
alpha(0.3)  %͸����
end    %��ѭ���������end

%% �Ӳ�ͬ�������ĸ����ڵ㣬�����������Ƕ��������ֵ
clc
clear all
fig=openfig('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\test\vdw gap-octahedrons_POS.fig')
pos=load('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\test\vdw gap-octahedrons_POS.txt'); %����txt����
y_min=min(pos(:,2));    %������������Ϊ���������߼�ֵ���֣���clear all���
range=5;    %����һ����Χ %ctrl+cǿ��ִ�У�������ֹ�����������
A=find(pos(:,2)<y_min+range);  %find�ҵ�y��Сֵ��Ӧ������
x_min=min(pos(A,1));    %�ҵ�yֵ��С��Ӧ��xֵ��Сֵ
[row,column]=find(pos(:,1)==x_min);  %find�ҵ�ֵ��Ӧ���к���
y=pos(row,2);   %��Сxֵ��Ӧ��yֵ
pos_first=[x_min,y] %��һ�����λ�á�

KNN=4; %�ҵ�һ����
Idx=knnsearch(pos,pos_first,'K',KNN);
four_spots=pos(Idx,:); %�ĸ�������
for i=1:1:length(four_spots);
    A=four_spots(i,:)-pos_first  %�ĸ������һ����Ĳ�ֵ    
end    
B=[atan2(A(1),A(5));atan2(A(2),A(6));atan2(A(3),A(7));atan2(A(4),A(8))] %�ҵ��ĸ���ĽǶ�
range_angle=0.2
find(B,pi/2-2<B<pi/2+0.2)









