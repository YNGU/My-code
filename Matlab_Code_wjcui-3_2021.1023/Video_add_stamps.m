%% add time stamps on in situ videos
clc
clear
%% 手动选择一张图加图表text，看效果
[f,p]=uigetfile('D:\++++++Data\Matlab Data\Exp_Data\GBT\Video add stamps\in situ TEM\key_frames\*.*','选择图像文件'); 
if f
A=imread(strcat(p,f));
A=imrotate(A,90);
imshow(A);              %展现图像
end
text(900,50,'0.5 s','horiz','center','color','w','fontsize',16,'FontName','Arial black')
% 3 是横坐标值，也是图像数组的列数值；4 是纵坐标值，也是数组的行数值
hold on 

%% read images 再批量设置导出带时间stamps的图片，再通过python的Filename_Batch operation.py改名字
folder='D:\++++++Data\Matlab Data\Exp_Data\GBT\Video add stamps\in situ TEM\key_frames\add_stamps\'; %文件夹后面要加/或者\
mfiles=dir([folder,'*.jpg']);  %dir输出文件夹中的数目
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
    ImageSum_R=imrotate(ImageSum_R,90);  %给图片转角度
    image_g_mat=fspecial('gaussian',[4 4],3);  %给一个高斯序列
    ImageSum_R=imfilter(ImageSum_R, image_g_mat); %运用高斯滤波
    ImageSum_R = insertText(ImageSum_R,[800 100],['t = ' num2str(0.5*i-0.5)  ' s'],...
        'FontSize',35,'Boxcolor',[255 255 255],'TextColor',[255 255 255],'Font','Arial',...
        'BoxOpacity',0);   %round(0.5*i-0.5)
    imagesc(ImageSum_R);  %看图
    axis image %
    png_file_name=['D:\++++++Data\Matlab Data\Exp_Data\GBT\Video add stamps\in situ TEM\key_frames\add_stamps\' name(1:end-3) 'jpg']; %输出位置
    imwrite(ImageSum_R,png_file_name);
end