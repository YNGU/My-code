%% clear
clear;
clc;
%% load mat file
load('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\Bi��Ge escape\30.mat');
%%

%% load position of txt files Fig.2
pos=load('D:\++++++Data\Matlab Data\Exp_Data\Test\16pos.txt'); %����txt����

folder='E:\'; %��һ���ļ��У�new_folder='D:\test'�������ļ���


f=figure; %����һ��fig�Ĵ���
f.PaperPosition=[1 1 7 14]; %Position��λ������[left��bottom��width��height]��ɣ�����������λ��
set(gca,'fontsize',9) %���������С
image=ReadDMFile('D:\++++++Data\Matlab Data\Exp_Data\Test\16pos.dm3'); %��DM�ļ�����ֻ������ߣ���ȥ��һ�伴��
image=double(image); 
imagesc(image); %������ image �е�������ʾΪһ��ͼ��
axis image off %���ֿ�߱Ȳ�ȡ��������
ylim([0 2048]); %��ǰ�������ͼ��� y ������
xlim([0 2048]);
colormap(Dan_Mumford); %��ɫ�ʱ��
hold all %�����Ļ�ͼ���������֮ǰ�Ļ�ͼ���������ס��ǰʹ�õ����ͺ���ɫ���Ӷ��µĻ�ͼ����ʹ�ò�ͬ����ɫ��������ʾ����
%plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k'); %pos(:,1)Ϊx,%MarkerEdgeColor���������ñ�ǵ�ı߿�������ɫ;MarkerFaceColor���������ñ�ǵ���ڲ����������ɫ;Markersize���������ñ�ǵ�Ĵ�С
%plot(pos(:,2)+1,pos(:,1)-1,'o','markersize',3,'markerfacecolor','w','markeredgecolor','none');�ƶ����Ӧ��λ��+1-1
x=pos(:,2); %��ȡ��ά����ĵ�һ�����ݣ�Ҳ����x
y=pos(:,1); %��ȡ��ά����ĵڶ������ݣ�Ҫ����ά�������������Ҫ����ֵ��:��ʾ
imagesc(image); %������ image �е�������ʾΪһ��ͼ��
axis image off %���ֿ�߱Ȳ�ȡ��������
ylim([700 1530]); %��ǰ�������ͼ��� y ������
xlim([500 720]);
colormap(Dan_Mumford); %��ɫ�ʱ��
hold all

KNN=6; % K���ڽ��㷨����ʵ�ַ������бȽϼ��׶���һ�ַ����㷨����6����

Idx=knnsearch(pos,pos,'K',KNN); %Idx = knnsearch(X,Y,Name,Value)��Idx = knnsearch(X,Y) ΪY�е�ÿ����ѯ�����X�е�����ھ�
Idx(:,1)=[]; %���þ����һ��Ϊ��ֵ��0��Ҳ����ȥ������Ϊ����ĵ� ��:,1��idx�ĵ�һ������

%Idx_angle=Idx-Idx; %��һ������ûɶ��
for i=1:1:length(pos)  %for��ѭ��������length����������Ԫ�ظ�����Ҳ����pos���ж��ٸ���
    for j=1:1:KNN-1  %���ڵ�������������Լ�
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
            '-r','linewidth',0.7);  %plot����
    end
    
end
%plot(x,y,'-c','linewidth',0.5) 
plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k');  %����
print(f,'-dpng','-r300',['D:\++++++Data\Matlab Data\Exp_Data\BT-Mn\Sb2Te3-Mn\Calatom\1-5.png']); 
%%  --
%pos=load('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\Bi��Ge escape\29pos.txt');
x=pos(:,2); %��ȡ��ά����ĵ�һ�����ݣ�Ҳ����x
y=pos(:,1); %��ȡ��ά����ĵڶ������ݣ�Ҫ����ά�������������Ҫ����ֵ��:��ʾ

%% save to mat file
save('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\Bi��Ge escape\30.mat');

%% �ԽǶ��ҵ�һ����
imagesc(image); %������ image �е�������ʾΪһ��ͼ��
axis image off %���ֿ�߱Ȳ�ȡ��������
ylim([0 800]); %��ǰ�������ͼ��� y ������
xlim([400 600]);
colormap(inferno); %��ɫ�ʱ��
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

%% �Ƕȵķֲ�ֱ��ͼ
index_ang=find(Idx_angle<150&Idx_angle>50) %�ҵ�С��170�������ж��ٸ�
ang=Idx_angle(index_ang) 
%hist(Idx_angle(:),100); %�ҵ��Ƕȷֲ�ֱ��ͼ
hist(ang(:),100);

%% load position of txt files Fig.3
pos=load('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\001_surface\12-pos-2.txt'); %����txt����

f=figure; %����һ��fig�Ĵ���
f.PaperPosition=[1 1 7 14]; %Position��λ������[left��bottom��width��height]��ɣ�����������λ��
set(gca,'fontsize',9) %���������С
%image=ReadDMFile('D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\001_surface\1-inferno.dm3'); %��DM�ļ�
image=double(image); 
imagesc(image); %������ image �е�������ʾΪһ��ͼ��
axis image off %���ֿ�߱Ȳ�ȡ��������
ylim([0 1024]); %��ǰ�������ͼ��� y ������
xlim([300 500]);
colormap(inferno); %��ɫ�ʱ��
hold all %�����Ļ�ͼ���������֮ǰ�Ļ�ͼ���������ס��ǰʹ�õ����ͺ���ɫ���Ӷ��µĻ�ͼ����ʹ�ò�ͬ����ɫ��������ʾ����
%plot(pos(:,2),pos(:,1),'o','markersize',2.5,'markerfacecolor','w','markeredgecolor','k'); %pos(:,1)Ϊx,%MarkerEdgeColor���������ñ�ǵ�ı߿�������ɫ;MarkerFaceColor���������ñ�ǵ���ڲ����������ɫ;Markersize���������ñ�ǵ�Ĵ�С
%plot(pos(:,2)+1,pos(:,1)-1,'o','markersize',3,'markerfacecolor','w','markeredgecolor','none');�ƶ����Ӧ��λ��+1-1
x=pos(:,2); %��ȡ��ά����ĵ�һ�����ݣ�Ҳ����x
y=pos(:,1); %��ȡ��ά����ĵڶ������ݣ�Ҫ����ά�������������Ҫ����ֵ��:��ʾ
imagesc(image); %������ image �е�������ʾΪһ��ͼ��
axis image off %���ֿ�߱Ȳ�ȡ��������
ylim([0 1024]); %��ǰ�������ͼ��� y ������
xlim([300 500]);
colormap(inferno); %��ɫ�ʱ��
hold all

KNN=6; % K���ڽ��㷨����ʵ�ַ������бȽϼ��׶���һ�ַ����㷨
Idx=knnsearch(pos,pos,'K',KNN); %Idx = knnsearch(X,Y,Name,Value)��Idx = knnsearch(X,Y) ΪY�е�ÿ����ѯ�����X�е�����ھ�
Idx(:,1)=[]; %���þ����һ��Ϊ��ֵ��0��Ҳ����ȥ������Ϊ����ĵ�

Idx_angle=Idx-Idx; %
for i=1:1:length(pos)  %for��ѭ��������
    for j=1:1:KNN-1  %
        Idx_angle(i,j)=atan2(pos(Idx(i,j),1)-pos(i,1),pos(Idx(i,j),2)-pos(i,2)); %����i��j֮��ĽǶ�
    end
end

range=0.5; %������Χ���������̫�̻���̫����Ҫ�Ĵ˲���
for i=1:1:length(pos) 
    % find the first vertex
    k1=find(abs(Idx_angle(i,:)-pi/2)<range); %��180�ȵĵ�
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
            '-k','linewidth',1);
    end
    
end
%plot(x,y,'-c','linewidth',0.5) 
plot(pos(:,2),pos(:,1),'o','markersize',3,'markerfacecolor','g','markeredgecolor','k');
print(f,'-dpng','-r300',['D:\Matlab Data\Exp Data\GBT data\GBT-CalAtom\001_surface\12plus-pos.png']);
%% pos

