%%
clc
clear
%% transform image file to inferno
image1=imread('D:\++++++Data\Matlab Data\Exp_Data\test\1550 11.4 Mx HAADF-DF4 Nano STEM 20220727 Diffraction.tif');
image1(:,:,4)=[]; %��ͨ����Ҫʱ����
image1=rgb2gray(image1);    %%��ͨ����Ҫʱ����
image1=double(image1);%��ͨ����Ҫʱ����
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
%caxis([50 250]) %����ͼ��ɫ�ʶԱȶ�
%% color_map������colormap��ȡΪpal�ļ�
clc;clear;
mycolorpoint=[[0 0 16];...    %6�����ص���Ϊ��ɫͼ�Ĳ�ֵ��������Բ�ֵ
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
load flujet  % matlab�Դ���flujetͼ
image(X);   %image,imagesc��imshow������չ��ͼ��
colormap(Dan_Mumford)
save('Dan_Mumford')
%%  
cmap2pal(Dan_Mumford) %ת��colormapΪ�������ļ�������origin�е�ɫ��
%%
cmap2pal(inferno)
%% Dan_Mumford ��ɫ ����֮ǰ������ǰ���save Dan��һ��
image1=imread('D:\++++++Data\Matlab Data\Exp_Data\BST-Mn\+Chapter2 Add\Bi2Te3-Mn-250 deg C 2h\eds��tip\SI EDS-HAADF 2041 31.8 kx Nano-HAADF.tif');

image1(:,:,4)=[];         %
image1=rgb2gray(image1);   %
image1=double(image1);
image2=0+(double(image1)-1)*55/255;   %ǰ���Ǳ��׵�ǿ��0��������100/450��ֵ��
image(image2);
colormap(Dan_Mumford);
axis image;
axis off;

%% ת���ļ���������tifͼ��Ϊinferno
clc
clear
Input_path = 'D:\++++++Data\1-1.QSTEMwDWfactor\GBTsCc1.5\test\'; %�����ļ���·�� 
Output_path = 'D:\++++++Data\1-1.QSTEMwDWfactor\GBTsCc1.5\test\'; %����ļ���·��
namelist = dir(strcat(Input_path, '*.tif')); %��ȡ���е�.tifͼƬ
len = length(namelist);
for i = 1:len
    name = namelist(i).name; %��ø�·���µ��ļ���
    image1=imread(strcat(Input_path, name));
    %image1(:,:,4)=[]; %��ͨ����Ҫʱ���� 
    image1=rgb2gray(image1);    %%��ͨ����Ҫʱ���� 
    image1=double(image1);%��ͨ����Ҫʱ���� 
    image2=40+(double(image1)-1)*155/255; 
    rgbImage = ind2rgb(image1, inferno(255)); 
    imshow(rgbImage);
    imagesc(rgbImage);
    colormap(inferno);
    axis image
    imwrite(rgbImage,[Output_path, name,'_inferno.tif']);  %ͼƬ����·����
end 

%% ת���ļ����Լ����ļ���������?