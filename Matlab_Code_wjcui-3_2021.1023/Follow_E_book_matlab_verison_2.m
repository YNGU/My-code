%% 2021.1211  图像的输入输出与jpg和tif的格式保存
image=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH01_Images\Fig0104.tif.');
imshow(image);
%figure, imshow(g);    %用figure保留第一张图，并显示第二张图
%imwrite(image,'D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH01_Images\1.jpg','quality',50) ...
    %存图像位jpg格式，加上‘quality’，后面跟上数值0-100表示图像质量，越小图像越粗糙
imwrite(image,'D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH01_Images\1.tif', 'compression', ...
    'none', 'resolution', [1024 1024])     %导出tif文件的标准格式，'compression', 'parameter', 'resolution', [rolres rowres]， parameter是com的压缩数量，colres和rowres是行与列的分辨率
%% 2021.1213   数组索引
v=[1 3 5 7 9];
v_2=v(2);   %()表示索引的位置,矩阵后面跟小括号表示索引
w=v.';    %转置矩阵
v(1:3);    %:的使用
v(1:end);  %end就是最后一个数值
v([1 3 5]);  %对应数组索引
v(1:2:end);    %设置步长
A=[1 2 3; 4 5 6; 7 8 9];
%A(2,:);   %:表示一行中的所有元素
A(2,1:end)  %与上一行等效果 
D=logical([1 0 0; 0 0 1; 0 0 0]);
A(D);
%% 函数句柄
f=@sin;   %@表面设置函数句柄 --第一类函数句柄-简单函数句柄
a=f(pi/4);
g=@(x) x.^2     % 第二个函数句柄-匿名函数句柄==有一个公式：@（input argument list） expression
g(2);
r=@(x,y) sqrt(x.^2 + y.^2);
r(3, 4);
%% 2021.1214 imadjust灰度转换
image=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH02_Images\Fig0203(a).tif');
image_contrast1=imadjust(image,[0 1],[1 0],1);    %imadjust函数中out的最大最小值相反即可灰度反转，值必须在0-1之间
image_contrast2=imcomplement(image);              %imcomplement也是灰度反转
figure(2);
subplot(2,2,1),imshow(image);
subplot(2,2,2),imshow(image_contrast1);
subplot(2,2,3),imshow(image_contrast2);
image2=imadjust(image,[0.5 0.75],[0 1],1);
subplot(2,2,4),imshow(image2);
image3=imadjust(image,[0 1],[0.25 0.75],1);
figure,imshow(image3);                         %要想打开另一个图片窗口，前面加figure，
%% 测试gamma值大于1小于1时的变化，默认值为1
image=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH02_Images\Fig0203(a).tif')
image_tran1=imadjust(image,[0 1],[0.2 0.9],1)
image_tran2=imadjust(image,[0 1],[0.2 0.9],2)
image_tran3=imadjust(image,[0 1],[0.2 0.9],3)
image_tran4=imadjust(image,[0 1],[0.2 0.9],0.7)
image_tran5=imadjust(image,[0 1],[0.2 0.9],0.3)
image_tran6=imadjust(image,[0 1],[0.2 0.9],0)      %gamma为0时就没有图像了
figure(1)
subplot(3,2,1),imshow(image_tran1);
subplot(3,2,2),imshow(image_tran2);
subplot(3,2,3),imshow(image_tran3);
subplot(3,2,4),imshow(image_tran4);
subplot(3,2,5),imshow(image_tran5);
subplot(3,2,6),imshow(image_tran6);
%% stretchlim(f)构成一个两元素向量，由一个低限和高限组成
%% log对比与对比度拉伸变换，改善一个高值区掩盖了低值的信号，为了将低值得信号也体现出来
f=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH02_Images\Fig0205(a).tif')
%imshow(f)
f1=mat2gray(f);                         %mat2gray将图片的pixel值限定在[0,1]之间
f2=im2uint8(f1);                           %im2uint8将图片变为unit8格式并且，pixel值限定在[0,255]之间
g=im2uint8(mat2gray(log(1+double(f))));      %log=loge；log2=log2
figure;
subplot(2,1,1),imshow(f);
subplot(2,1,2),imshow(g);
%% 小数可称为浮点数，作为一个小数3.14。如果使用指数表现形式的话（3.14E0）小数点可以任意浮动
%% 2021.1216 图像直方图+++加上分号之后就不会在下面命令行窗口运行，直接计算省了时间，
f=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH02_Images\Fig0205(a).tif'); 
h=imhist(f,25);     %imhist(f,N) N是把整个函数分为几点，类似下面的N
horz=linspace(0,255,25);    %linspace(x1,x2,N)中x1、x2、N分别为起始值、中止值、元素的个数。N默认点数为100。
bar(horz,h);       %bar,stem,plot是三种其他的画法
stem(horz,h);
plot(horz,h);
%imhist(f);   %imhist直接可以展示图片的灰度分布
%figure(1);
%subplot(2,2,1),imhist(f);
%subplot(2,2,2),bar(horz,h);
%subplot(2,2,3),stem(horz,h);
%subplot(2,2,4),plot(horz,h);     %可以直接在subplot，后面放函数，不用imshow（）也可以。
%%  xy轴的范围，间隔，标记
imhist(f);
axis([0 100 0 60000]);               %axis设置xy轴的最小值和最大值
set(gca, 'xtick', 0:10:100);
set(gca, 'ytick', 0:20000:60000);      % 设置x，y轴的间隔大小
xlabel('x axis', 'fontsize', 20);
ylabel('y axis', 'fontsize', 20);
title('gray distubrition','fontsize',20);   %设置题目
ylim('auto');
xlim('auto');

%%

