%% add time stamps on in situ videos
clc
clear
%% �ֶ�ѡ��һ��ͼ��ͼ��text����Ч��
[f,p]=uigetfile('D:\++++++Data\Matlab Data\Exp_Data\GBT\Video add stamps\in situ TEM\key_frames\*.*','ѡ��ͼ���ļ�'); 
if f
A=imread(strcat(p,f));
A=imrotate(A,90);
imshow(A);              %չ��ͼ��
end
text(900,50,'0.5 s','horiz','center','color','w','fontsize',16,'FontName','Arial black')
% 3 �Ǻ�����ֵ��Ҳ��ͼ�����������ֵ��4 ��������ֵ��Ҳ�����������ֵ
hold on 

%% read images ���������õ�����ʱ��stamps��ͼƬ����ͨ��python��Filename_Batch operation.py������
folder='D:\++++++Data\Matlab Data\Exp_Data\GBT\Video add stamps\in situ TEM\key_frames\add_stamps\'; %�ļ��к���Ҫ��/����\
mfiles=dir([folder,'*.jpg']);  %dir����ļ����е���Ŀ
l=length(mfiles);    
average_intensity=zeros(l,1); 
max_intensity=zeros(l,1);
min_intensity=zeros(l,1);
%c_high=0.15
for i=1:1:l  
    lname=length(mfiles(i).name);
    name=mfiles(i).name(1:lname);
    jpg_file_name=[folder name];
    ImageSum_R=im2double(imread(jpg_file_name));
    ImageSum_R=imrotate(ImageSum_R,90);  %��ͼƬת�Ƕ�
    image_g_mat=fspecial('gaussian',[4 4],3);  %��һ����˹����
    ImageSum_R=imfilter(ImageSum_R, image_g_mat); %���ø�˹�˲�
    ImageSum_R = insertText(ImageSum_R,[800 100],['t = ' num2str(0.5*i-0.5)  ' s'],...
        'FontSize',35,'Boxcolor',[255 255 255],'TextColor',[255 255 255],'Font','Arial',...
        'BoxOpacity',0);   %round(0.5*i-0.5)
    imagesc(ImageSum_R);  %��ͼ
    axis image %
    png_file_name=['D:\++++++Data\Matlab Data\Exp_Data\GBT\Video add stamps\in situ TEM\key_frames\add_stamps\' name(1:end-3) 'jpg']; %���λ��
    imwrite(ImageSum_R,png_file_name);
end