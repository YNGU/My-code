%% ZnO-Al2O3_Xupengyu_2022.0214 ʸ����ϵ
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
data_per=xlsread('D:\++++++Data\1\per.xlsx','sheet1','B1:B770');  %xlsread��excel������,excel���ڶ�ȡ·��
data_def=xlsread('D:\++++++Data\1\vsb.xlsx','sheet1','B1:B770');  
data=data_per-data_def
%% a�汾�ȶ��Բ��ߣ�b���ȶ��ԱȽϸ�---appdesigner��app���
a=1+1;   
b=1.*10;    %�����ĳ˳�����.����,����ƽ����,�η�����Ҫ��
c=b./a
%% �����matlabϵ�н̳� P1-P2____number
%command window
 %help+function / ���Ͻ�search����Ҫ�ĺ���
cos(((1+2+3+4)^3/5)^0.5);
sin(sqrt(pi))+log(tan(1));
2^(3.5*1.7);
exp(sin(10));
cos(pi):
sin(ans):
  %Variables����
  A=10; 
  a=20;
  A2=10;   %variable/data type, who��������Щ����,whos�����еı�����types/
  %ans,i,j complex number,inf,eps 2.2204e-016,NaN not a
  %number,pi,���Ѿ��к����variables,��iskeywords
  cos='This string.'; %�����cos�����variables,����cos������
  cos(8);  %�����cos(8)��ָcos�еڰ˸���ĸ;variables��built-in function�����ȼ�����
 %format long / shortE/longE/bank/hex/rat
%% �����matlabϵ�н̳� P1-P2____array
a=[1 2 3 4]; %������ row vector
b=[1; 2; 3; 4]; %������ column vector
c=a*b;
d=b*a;
A=[1, 21, 6; 5, 17, 9; 31, 2, 7];
a(3);
b(4);  %array indexing �Ҿ������������
A(1, 2);   %�Ҿ����λ����ֵ
A(2);
A(8);
A([1 3 5]);
A([1 3; 1 3]);
A([1 3], [1 3]); %ǰ��ķ��������ҵ�һ��������У�����ĵ�һ�к͵����С�
A(1,2)=16;
A([1,2],[2,3])=[0 0; 0 0]
A(3,:)={}; %��һ��ȫ��Ϊ0
F=[A B]; %�������
A=[1 2 3; 4 5 4; 9 8 7]
B=[3 3 3; 2 4 9; 1 3 1]   %+��-��*��/  ����ļӼ��˳��㷨,.*����Ƕ�Ӧλ�����,A/B=A*inv(B);A./B=��Ӧ��λ����ֵ�����'ת�þ���
a=2
A^a
A.^a
%some special Matrix
linspace(0 13 5)
eye(3)
zeros(4,2)
ones(4,2)
diag([2 3 4]) %�Խ��߾���
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
pos=load('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\Te1-Sb-Te1_octahedrons\Te1-Sb-Te1_octahedrons-Pos-28.txt'); %����txt����

folder='E:\'; %��һ���ļ��У�new_folder='D:\test'�������ļ���


f=figure; %����һ��fig�Ĵ���
f.PaperPosition=[1 1 7 14]; %Position��λ������[left��bottom��width��height]��ɣ�����������λ��
set(gca,'fontsize',9) %���������С
image=ReadDMFile('D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\1607 2.9 Mx STEM Diffraction HAADF 20211107 Nano-G.dm3'); %��DM�ļ�����ֻ������ߣ���ȥ��һ�伴��
image=double(image); 
%imagesc(image); %������ image �е�������ʾΪһ��ͼ��
%axis image off %���ֿ�߱Ȳ�ȡ��������
%ylim([0 2048]); %��ǰ�������ͼ��� y ������
%xlim([0 2048]);
%colormap(inferno); %��ɫ�ʱ��
%hold all %�����Ļ�ͼ���������֮ǰ�Ļ�ͼ���������ס��ǰʹ�õ����ͺ���ɫ���Ӷ��µĻ�ͼ����ʹ�ò�ͬ����ɫ��������ʾ����
%plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k'); %pos(:,1)Ϊx,%MarkerEdgeColor���������ñ�ǵ�ı߿�������ɫ;MarkerFaceColor���������ñ�ǵ���ڲ����������ɫ;Markersize���������ñ�ǵ�Ĵ�С
%plot(pos(:,2)+1,pos(:,1)-1,'o','markersize',3,'markerfacecolor','w','markeredgecolor','none');�ƶ����Ӧ��λ��+1-1
x=pos(:,2); %��ȡ��ά����ĵ�һ�����ݣ�Ҳ����x
y=pos(:,1); %��ȡ��ά����ĵڶ������ݣ�Ҫ����ά�������������Ҫ����ֵ��:��ʾ
imagesc(image); %������ image �е�������ʾΪһ��ͼ��
axis image off %���ֿ�߱Ȳ�ȡ��������
ylim([0 2048]); %��ǰ�������ͼ��� y ������
xlim([0 2048]);
colormap(inferno); %��ɫ�ʱ��
hold all

KNN=6; % K���ڽ��㷨����ʵ�ַ������бȽϼ��׶���һ�ַ����㷨����6����

Idx=knnsearch(pos,pos,'K',KNN); %Idx = knnsearch(X,Y,Name,Value)��Idx = knnsearch(X,Y) ΪY�е�ÿ����ѯ�����X�е�����ھ�
Idx(:,1)=[]; %���þ����һ��Ϊ��ֵ��0��Ҳ����ȥ������Ϊ����ĵ� ��:,1��idx�ĵ�һ������

%Idx_angle=Idx-Idx; %��һ������ûɶ��
for i=1:1:length(pos)  %for��ѭ��������length����������Ԫ�ظ������趨i��pos�еĵڼ���λ��
    for j=1:1:KNN-1  %���ڵ�������������Լ����趨jΪKNN�еĵ�
        Idx_angle(i,j)=atan2(pos(Idx(i,j),1)-pos(i,1),pos(Idx(i,j),2)-pos(i,2)); %����i��j֮��ĽǶ�,atan2(x,y)��������
    end
end

range=0.5; %������Χ���������̫�̻���̫����Ҫ�Ĵ˲���
for i=1:1:length(pos) 
    % find the first vertex
    k1=find(abs(Idx_angle(i,:)-pi/2)<range); %��180�ȵĵ�, abs����һ�����ľ���ֵ
    if isempty(k1) %�ж��Ƿ�Ϊ��
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
            '-r','linewidth',0.7);  %plot���ߣ����ҳ��İ������ĸ�������
    end
    if mod(i,1)==0
        c=rand(28); %�����ɫ����
        patch([pos(i,2) pos(index1,2) pos(index2,2) pos(index3,2) pos(i,2)],...
            [pos(i,1) pos(index1,1) pos(index2,1) pos(index3,1) pos(i,1)],...
            c);  %patch����ĸ����ƽ���ı���Ϊ��ɫ
        colormap(inferno)
    end
    
end
plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k');  %��ԭ��λ�û���
print(f,'-dpng','-r300',['D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\1-5-test.png']); 
%% ����ĸ����������
patch(pos(:,1),pos(:,2),'red');
imagesc(imgae)
hold all;


%% forѭ�����
clear;
sum=0;
for i=1:100;
    sum=sum+i;
end
sum=sum %�򵥵�forѭ�����

clear;
sum=0;
for i=1:2:100;   %2�����õĲ���
    sum=sum+i;
end
sum=sum    %���ϲ�����forѭ�����

clear;
A=rand(1,4);
for i=A
    i
end            %�������ı�����ָ����ĳ������·�ߣ����ζ�����ÿ��������һ���ҽ���һ�η���

clear;
for i=1:1:5
    for j=1:1:3
        A=i*j
    end
end

%% if ѭ�����
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

 
  
  