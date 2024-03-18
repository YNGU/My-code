%%
clc
clear
%% transform image file to inferno
image1=imread('D:\++++++Data\Matlab Data\Exp_Data\test\1550 11.4 Mx HAADF-DF4 Nano STEM 20220727 Diffraction.tif');
image1(:,:,4)=[]; %四通道必要时再用
image1=rgb2gray(image1);    %%四通道必要时再用
image1=double(image1);%四通道必要时再用
image2=40+(double(image1)-1)*50/255;
rgbImage = ind2rgb(image1, inferno(255));  
imshow(rgbImage);
imagesc(rgbImage);
colormap(inferno);
axis image
imwrite(rgbImage,'D:\++++++Data\Matlab Data\Exp_Data\test\1550 11.4 Mx HAADF-DF4 Nano STEM 20220727 Diffraction_inferno.tif');
%% transform image file to hot
image1=imread('D:\++++++Data\Matlab Data\Exp_Data\test\1550 11.4 Mx HAADF-DF4 Nano STEM 20220727 Diffraction.tif');
image1(:,:,4)=[];
image1=rgb2gray(image1);
image1=double(image1);
image2=40+(double(image1)-1)*155/255;
rgbImage = ind2rgb(image1, hot(255));   %Dan_Mumford/hot
imshow(rgbImage);
imagesc(rgbImage);
colormap(hot);
axis image;
imwrite(rgbImage,'D:\++++++Data\Matlab Data\Exp_Data\test\1550 11.4 Mx HAADF-DF4 Nano STEM 20220727 Diffraction_hot.tif');
%%
%imagesc(image1)
%axis image
%colormap(inferno)
%caxis([50 250]) %控制图像色彩对比度
%% color_map制作与colormap提取为pal文件
clc;clear;
mycolorpoint=[[0 0 16];...    %6个像素点作为颜色图的插值点进行线性插值
    [8 69 99];...
    [57 174 156];...
    [198 243 99];...
    [222 251 123];...
    [239 255 190]]
mycolorposition=[1 11 33 50 57 64];
mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:64,'linear','extrap');
mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:64,'linear','extrap');
mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:64,'linear','extrap');
Dan_Mumford=[mycolormap_r',mycolormap_g',mycolormap_b']/255;
Dan_Mumford=round(Dan_Mumford*10^4)/10^4
load flujet  % matlab自带的flujet图
image(X);   %image,imagesc和imshow都可以展现图像
colormap(Dan_Mumford)
save('Dan_Mumford')
%%  
cmap2pal(Dan_Mumford) %转化colormap为二进制文件，用于origin中调色板
%%
cmap2pal(inferno)
%% Dan_Mumford 颜色 运行之前先运行前面的save Dan那一步
image1=imread('D:\++++++Data\Matlab Data\Exp_Data\BST-Mn\+Chapter2 Add\Bi2Te3-Mn-250 deg C 2h\eds—tip\SI EDS-HAADF 2041 31.8 kx Nano-HAADF.tif');

image1(:,:,4)=[];         %
image1=rgb2gray(image1);   %
image1=double(image1);
image2=0+(double(image1)-1)*55/255;   %前面是背底的强度0，变后面的100/450的值，
image(image2);
colormap(Dan_Mumford);
axis image;
axis off;

%% 转变文件夹中所有tif图像为inferno
clc
clear
Input_path = 'D:\++++++Data\1-1.QSTEMwDWfactor\GBTsCc1.5\test\'; %输入文件夹路径 
Output_path = 'D:\++++++Data\1-1.QSTEMwDWfactor\GBTsCc1.5\test\'; %输出文件夹路径
namelist = dir(strcat(Input_path, '*.tif')); %获取所有的.tif图片
len = length(namelist);
for i = 1:len
    name = namelist(i).name; %获得该路径下的文件名
    image1=imread(strcat(Input_path, name));
    %image1(:,:,4)=[]; %四通道必要时再用 
    image1=rgb2gray(image1);    %%四通道必要时再用 
    image1=double(image1);%四通道必要时再用 
    image2=40+(double(image1)-1)*155/255; 
    rgbImage = ind2rgb(image1, inferno(255)); 
    imshow(rgbImage);
    imagesc(rgbImage);
    colormap(inferno);
    axis image
    imwrite(rgbImage,[Output_path, name,'_inferno.tif']);  %图片存续路径名
end 

%% 转变文件夹以及子文件夹中所有?