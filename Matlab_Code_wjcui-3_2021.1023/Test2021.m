%% 2021.0829 读取与保存图片
image=imread('D:\Matlab Data\Exp Data\GBT data\Trans inferno imagesSTEM HAADF 1333 20201104 31.8 kx-inferno.tif');
%2048x2018x3前面2048先2048是像素，x3说明是彩色图片代表RGB三个颜色通道
imshow(image); 
image_R=image(:,:,1); %找到图像中的R通道
image_G=image(:,:,2);
image_B=image(:,:,3);

figure(2);  %创建figure2显示一下图像的对话框
subplot(2,2,1),imshow(image_R);  %subplot将对话框分成2*2个子部分，1并在第一部分显示图像
%subplot是将多个图画到一个平面上的工具，subplot(m,n,p)生成m*n个子图，当前激活第p个子图
subplot(2,2,2),imshow(image_G);
subplot(2,2,3),imshow(image_B);
%subplot(2,2,4),imshow(image);

image_gray=rgb2gray(image) %RGB图像转变为灰度图
imwrite(image_gray,'D:\Matlab Data\Code Sum\Matlab_Code_wjcui-3\Test2021\1.png')  
%保存图像为png格式
%% 拉伸/压缩图像
image_resize=imresize(image,[200,200]); %imresize改变图像像素的数量，拉伸压缩
image_rotate=imrotate(image,30); %对图像进行旋转操作，数字表示旋转的度数
figure(2);
imshow(image_rotate); 
subplot(2,1,1),imshow(image);
subplot(2,1,2),imshow(image_rotate);
image_flip=flipdim(image_resize,1)  %对图像进行镜像，1为左右，2为上下对称
imshow(image_flip)
%% 图像滤波——滤波就是去除图像中的噪声
%对于数字图像信号，噪声表为或大或小的极值，这些极值通过加减作用于图像像素的真实灰度值上，
%对图像造成亮、暗点干扰，极大降低了图像质量，影响图像复原、分割、特征提取、图像识别等后继工作的进行。
%高斯噪声：这种噪声符合正态分布
%椒盐噪声，均匀噪声，伽马噪声等等
%算子（Operator）：对一个或几个算数进行计算处理，如加、减、乘、除等。为了方便，我们也可以把函数当作算子
%Sobel_Operator=(1 0 -1; 2 0 -2; 1 0 -1),  Sobel算子：用来检测图像边缘的滤波模板

%中值滤波：中值滤波由Turky在1971年提出，最初用于时间序列分析，后来被用于图像处理；基本原理是把图像或序列中心点位置的值用该域的中值替代，曾被认为是非线性滤波的代表
%在原始数据的周围选取一定的区域，将这个区域内所有的数都排序，选取排序后中间的一个数来代替原始数据。
%在采用中值滤波的时候一定要使用奇数个像素作为模板，不然中间的数就不会是一个数，是的滤波时候出现问题

%高斯滤波：
%% 影响科研的只有两门学科：数学和英语。英语决定你在某个领域的广度；数学觉得你在某个领域的深度。
%% 常用图像噪声引入
image=imread('D:\Matlab Data\Code Sum\Matlab_Code_wjcui-3\Test2021\origin.jpg')
image_gray=rgb2gray(image);  %RGB转灰度
image_gaussian=imnoise(image_gray, 'gaussian', 0, 39*39/(255*255));  %给图像加高斯噪声，后面两个数是高斯噪声（函数）的均值和σ
image_salt_pepper=imnoise(image_gray, 'salt & pepper', 0.05);  %加椒盐噪声，后面的参数是椒盐噪声的密度
subplot(1,2,1),imshow(image_gaussian);
subplot(1,2,2),imshow(image_salt_pepper);

%% 编程均值滤波-椒盐滤波：主要分为4个步骤：生成滤波模板、延拓被滤波图像、像素位置遍历滤波、裁剪滤波后图像
n=5;  %滤波模板大小
[height, width]=size(image_salt_pepper);   %输入图像是p×q的,且p>n,q>n  
image_salt_pepper2=zeros(height+n-1,width+n-1); %存放延拓后的图片
for i=1+(n-1)/2:1+(n-1)/2+height-1
    for j=1+(n-1)/2:1+(n-1)/2+width-1
        image_salt_pepper2(i,j)=image_salt_pepper(i-(n-1)/2,j-(n-1)/2);
    end
end
image_salt_pepper3=image_salt_pepper2;   %灰度值更改后会对后续的滤波产生影响，因此需要一个不变的图片
a(1:n,1:n)=1;
for i=1:height  
   for j=1:width  
       c=image_salt_pepper3(i:i+(n-1),j:j+(n-1)).*a; %取出x1中从(i,j)开始的n行n列元素与模板相乘  
       s=sum(sum(c));                 %求c矩阵中各元素之和  
       image_salt_pepper2(i+(n-1)/2,j+(n-1)/2)=s/(n*n); %将与模板运算后的各元素的均值赋给模板中心位置的元素  
   end  
end  
image_salt_pepper2=uint8(image_salt_pepper2((n-1)/2+1:height+(n-1)/2,(n-1)/2+1:width+(n-1)/2));%通过计算，此处只是一个矩阵，用uint8将矩阵变成图像
imshow(image_salt_pepper2);
title('wechat均值滤波')

%% 编程中值滤波：分为3个步骤：延拓被滤波图像、像素位置遍历滤波、裁剪滤波后图像。不需要生成滤波模板
for i=1:height  
     for j=1:width 
         c=image_salt_pepper3(i:i+(n-1),j:j+(n-1)); %取出x1中从(i,j)开始的n行n列元素,即模板(n×n的)  
         e=c(1,:);      %是c矩阵的第一行  
         for u=2:n  
            e=[e,c(u,:)];     %将c矩阵变为一个行矩阵      
         end  
         mm=median(e);      %mm是中值  
        image_salt_pepper2(i+(n-1)/2,j+(n-1)/2)=mm;   %将模板各元素的中值赋给模板中心位置的元素  
     end
end
image_salt_pepper2=uint8(image_salt_pepper2((n-1)/2+1:height+(n-1)/2,(n-1)/2+1:width+(n-1)/2));
imshow(image_salt_pepper2);
title('wechat中值滤波')
%% 生成滤波模板：先生成滤波模板，再进行滤波函数的使用
 image_g_mat=fspecial('gaussian',[5 5],3); %生成高斯序列
 % image_g_mat=[0.003 0.013 0.022 0.013 0.003;
 %      0.013 0.059 0.097 0.059 0.013;
 %      0.022 0.097 0.159 0.097 0.022;
 %      0.013 0.059 0.097 0.059 0.013;
 %      0.003 0.013 0.022 0.013 0.003];
image_avr_mat=fspecial('average', 5);  %生成均值滤波序列
%A_avr_mat=[0.1111 0.1111 0.1111 0.1111;
%            0.1111 0.1111 0.1111 0.1111;
%            0.1111 0.1111 0.1111 0.1111;
%            0.1111 0.1111 0.1111 0.1111;]
%处理高斯噪声
image_g=imfilter(image_gaussian, image_g_mat);   %用生成的高斯序列进行滤波
image_m=medfilt2(image_gaussian,[9 9]);   %中值滤波，后面的参数为滤波模板的尺寸
image_a=imfilter(image_gaussian, image_avr_mat);
subplot(2,2,1),imshow(image_gaussian),title('含有高斯噪声的图像');
subplot(2,2,2),imshow(image_g),title('高斯滤波处理高斯噪声后的图像');
subplot(2,2,3),imshow(image_m),title('中值滤波处理高斯噪声后的图像');
subplot(2,2,4),imshow(image_a),title('均值滤波处理高斯噪声后的图像');
%无论采用哪种滤波方式，都无法完全的去除掉图像中的噪声，只能最大程度上的减小这种噪声对图像的影响。
%三种噪声各有自己的长出，而且每种滤波方式都可以通过修改滤波模板来调整滤波的效果，在实际的应用中需要根据情况灵活的选择，并加以运用
%%
clc
clear
%%
A=magic(6) 
imshow(A,'InitialMagnification','fit') %以适合窗口大小显示A图

%%  寻找最近邻点
KNN=6;
[Idx,d]=knnsearch(pos,pos,'K',KNN);  %找到两个pos矩阵的相邻的点idx，6个，并d为它们对应的距离
Idx=knnsearch(pos,pos,'K',KNN);  %找到相邻的点idx，不记录记录
Idx = knnsearch(X,Y)  %找y中与x最近邻点，一个存于idx中
%% 手动输入a，c角度和l1-l4
ang_a=[];
ang_c=[];
dis_l1=[];
dis_l2=[];
dis_l3=[];
dis_l4=[];
%%  求四边形面积1743
S1=0.5.*sind(ang_a).*(0.024.*dis_l1).*(0.024.*dis_l4)  %单位是nm2,,矩阵数对应相乘用点乘(.*)
S2=0.5.*sind(ang_c).*(0.024.*dis_l2).*(0.024.*dis_l3)   %pxiel=0.024*0.024
S=S1+S2 
%% 求四边形面积1754
S1=0.5.*sind(ang_a).*(0.017.*dis_l1).*(0.017.*dis_l4)  %单位是nm2,,矩阵数对应相乘用点乘(.*)
S2=0.5.*sind(ang_c).*(0.017.*dis_l2).*(0.017.*dis_l3)   %pxiel=0.017*0.017
S=S1+S2 

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
%% y=k*x+b
x1=-79.29
y1=-0.306
x2=-27.60
y2=-0.9751
k=(y1+y2)/(x1+x2);
b=k*x1-y1;
y_x0=b
x_y0=-b/k

%% 
dis=load('C:\Users\Lenovo\Desktop\dis.txt')

%% cal
M=27*1.6606e-27  % Al的原子质量=相对原子质量*1.6606e-27 kg
E=80000    %单位keV
m=9.1e-31   %静止的单个电子质量，单位kg
c=3e8    %光速m/s
E1=2*M*E*(E+2*m*(c^2))
E2=((M+m)^2)*(c^2)+2*M*E
Em=E1/E2

%% 
M=72.6*1.6606e-27 %原子质量
E=80000*1.602e-19 %J
m=9.1e-31   %静止的单个电子质量，单位kg
c=300000000    %光速m/s

Em=((2*E*(E+2*m*(c^2)))/(M*(c^2)))/1.602e-19



