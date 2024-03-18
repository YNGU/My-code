%% 2021.0829 ��ȡ�뱣��ͼƬ
image=imread('D:\Matlab Data\Exp Data\GBT data\Trans inferno imagesSTEM HAADF 1333 20201104 31.8 kx-inferno.tif');
%2048x2018x3ǰ��2048��2048�����أ�x3˵���ǲ�ɫͼƬ����RGB������ɫͨ��
imshow(image); 
image_R=image(:,:,1); %�ҵ�ͼ���е�Rͨ��
image_G=image(:,:,2);
image_B=image(:,:,3);

figure(2);  %����figure2��ʾһ��ͼ��ĶԻ���
subplot(2,2,1),imshow(image_R);  %subplot���Ի���ֳ�2*2���Ӳ��֣�1���ڵ�һ������ʾͼ��
%subplot�ǽ����ͼ����һ��ƽ���ϵĹ��ߣ�subplot(m,n,p)����m*n����ͼ����ǰ�����p����ͼ
subplot(2,2,2),imshow(image_G);
subplot(2,2,3),imshow(image_B);
%subplot(2,2,4),imshow(image);

image_gray=rgb2gray(image) %RGBͼ��ת��Ϊ�Ҷ�ͼ
imwrite(image_gray,'D:\Matlab Data\Code Sum\Matlab_Code_wjcui-3\Test2021\1.png')  
%����ͼ��Ϊpng��ʽ
%% ����/ѹ��ͼ��
image_resize=imresize(image,[200,200]); %imresize�ı�ͼ�����ص�����������ѹ��
image_rotate=imrotate(image,30); %��ͼ�������ת���������ֱ�ʾ��ת�Ķ���
figure(2);
imshow(image_rotate); 
subplot(2,1,1),imshow(image);
subplot(2,1,2),imshow(image_rotate);
image_flip=flipdim(image_resize,1)  %��ͼ����о���1Ϊ���ң�2Ϊ���¶Գ�
imshow(image_flip)
%% ͼ���˲������˲�����ȥ��ͼ���е�����
%��������ͼ���źţ�������Ϊ����С�ļ�ֵ����Щ��ֵͨ���Ӽ�������ͼ�����ص���ʵ�Ҷ�ֵ�ϣ�
%��ͼ���������������ţ����󽵵���ͼ��������Ӱ��ͼ��ԭ���ָ������ȡ��ͼ��ʶ��Ⱥ�̹����Ľ��С�
%��˹��������������������̬�ֲ�
%��������������������٤�������ȵ�
%���ӣ�Operator������һ���򼸸��������м��㴦����ӡ������ˡ����ȡ�Ϊ�˷��㣬����Ҳ���԰Ѻ�����������
%Sobel_Operator=(1 0 -1; 2 0 -2; 1 0 -1),  Sobel���ӣ��������ͼ���Ե���˲�ģ��

%��ֵ�˲�����ֵ�˲���Turky��1971��������������ʱ�����з���������������ͼ��������ԭ���ǰ�ͼ����������ĵ�λ�õ�ֵ�ø������ֵ�����������Ϊ�Ƿ������˲��Ĵ���
%��ԭʼ���ݵ���Χѡȡһ�������򣬽�������������е���������ѡȡ������м��һ����������ԭʼ���ݡ�
%�ڲ�����ֵ�˲���ʱ��һ��Ҫʹ��������������Ϊģ�壬��Ȼ�м�����Ͳ�����һ�������ǵ��˲�ʱ���������

%��˹�˲���
%% Ӱ����е�ֻ������ѧ�ƣ���ѧ��Ӣ�Ӣ���������ĳ������Ĺ�ȣ���ѧ��������ĳ���������ȡ�
%% ����ͼ����������
image=imread('D:\Matlab Data\Code Sum\Matlab_Code_wjcui-3\Test2021\origin.jpg')
image_gray=rgb2gray(image);  %RGBת�Ҷ�
image_gaussian=imnoise(image_gray, 'gaussian', 0, 39*39/(255*255));  %��ͼ��Ӹ�˹�����������������Ǹ�˹�������������ľ�ֵ�ͦ�
image_salt_pepper=imnoise(image_gray, 'salt & pepper', 0.05);  %�ӽ�������������Ĳ����ǽ����������ܶ�
subplot(1,2,1),imshow(image_gaussian);
subplot(1,2,2),imshow(image_salt_pepper);

%% ��̾�ֵ�˲�-�����˲�����Ҫ��Ϊ4�����裺�����˲�ģ�塢���ر��˲�ͼ������λ�ñ����˲����ü��˲���ͼ��
n=5;  %�˲�ģ���С
[height, width]=size(image_salt_pepper);   %����ͼ����p��q��,��p>n,q>n  
image_salt_pepper2=zeros(height+n-1,width+n-1); %������غ��ͼƬ
for i=1+(n-1)/2:1+(n-1)/2+height-1
    for j=1+(n-1)/2:1+(n-1)/2+width-1
        image_salt_pepper2(i,j)=image_salt_pepper(i-(n-1)/2,j-(n-1)/2);
    end
end
image_salt_pepper3=image_salt_pepper2;   %�Ҷ�ֵ���ĺ��Ժ������˲�����Ӱ�죬�����Ҫһ�������ͼƬ
a(1:n,1:n)=1;
for i=1:height  
   for j=1:width  
       c=image_salt_pepper3(i:i+(n-1),j:j+(n-1)).*a; %ȡ��x1�д�(i,j)��ʼ��n��n��Ԫ����ģ�����  
       s=sum(sum(c));                 %��c�����и�Ԫ��֮��  
       image_salt_pepper2(i+(n-1)/2,j+(n-1)/2)=s/(n*n); %����ģ�������ĸ�Ԫ�صľ�ֵ����ģ������λ�õ�Ԫ��  
   end  
end  
image_salt_pepper2=uint8(image_salt_pepper2((n-1)/2+1:height+(n-1)/2,(n-1)/2+1:width+(n-1)/2));%ͨ�����㣬�˴�ֻ��һ��������uint8��������ͼ��
imshow(image_salt_pepper2);
title('wechat��ֵ�˲�')

%% �����ֵ�˲�����Ϊ3�����裺���ر��˲�ͼ������λ�ñ����˲����ü��˲���ͼ�񡣲���Ҫ�����˲�ģ��
for i=1:height  
     for j=1:width 
         c=image_salt_pepper3(i:i+(n-1),j:j+(n-1)); %ȡ��x1�д�(i,j)��ʼ��n��n��Ԫ��,��ģ��(n��n��)  
         e=c(1,:);      %��c����ĵ�һ��  
         for u=2:n  
            e=[e,c(u,:)];     %��c�����Ϊһ���о���      
         end  
         mm=median(e);      %mm����ֵ  
        image_salt_pepper2(i+(n-1)/2,j+(n-1)/2)=mm;   %��ģ���Ԫ�ص���ֵ����ģ������λ�õ�Ԫ��  
     end
end
image_salt_pepper2=uint8(image_salt_pepper2((n-1)/2+1:height+(n-1)/2,(n-1)/2+1:width+(n-1)/2));
imshow(image_salt_pepper2);
title('wechat��ֵ�˲�')
%% �����˲�ģ�壺�������˲�ģ�壬�ٽ����˲�������ʹ��
 image_g_mat=fspecial('gaussian',[5 5],3); %���ɸ�˹����
 % image_g_mat=[0.003 0.013 0.022 0.013 0.003;
 %      0.013 0.059 0.097 0.059 0.013;
 %      0.022 0.097 0.159 0.097 0.022;
 %      0.013 0.059 0.097 0.059 0.013;
 %      0.003 0.013 0.022 0.013 0.003];
image_avr_mat=fspecial('average', 5);  %���ɾ�ֵ�˲�����
%A_avr_mat=[0.1111 0.1111 0.1111 0.1111;
%            0.1111 0.1111 0.1111 0.1111;
%            0.1111 0.1111 0.1111 0.1111;
%            0.1111 0.1111 0.1111 0.1111;]
%�����˹����
image_g=imfilter(image_gaussian, image_g_mat);   %�����ɵĸ�˹���н����˲�
image_m=medfilt2(image_gaussian,[9 9]);   %��ֵ�˲�������Ĳ���Ϊ�˲�ģ��ĳߴ�
image_a=imfilter(image_gaussian, image_avr_mat);
subplot(2,2,1),imshow(image_gaussian),title('���и�˹������ͼ��');
subplot(2,2,2),imshow(image_g),title('��˹�˲������˹�������ͼ��');
subplot(2,2,3),imshow(image_m),title('��ֵ�˲������˹�������ͼ��');
subplot(2,2,4),imshow(image_a),title('��ֵ�˲������˹�������ͼ��');
%���۲��������˲���ʽ�����޷���ȫ��ȥ����ͼ���е�������ֻ�����̶��ϵļ�С����������ͼ���Ӱ�졣
%�������������Լ��ĳ���������ÿ���˲���ʽ������ͨ���޸��˲�ģ���������˲���Ч������ʵ�ʵ�Ӧ������Ҫ�����������ѡ�񣬲���������
%%
clc
clear
%%
A=magic(6) 
imshow(A,'InitialMagnification','fit') %���ʺϴ��ڴ�С��ʾAͼ

%%  Ѱ������ڵ�
KNN=6;
[Idx,d]=knnsearch(pos,pos,'K',KNN);  %�ҵ�����pos��������ڵĵ�idx��6������dΪ���Ƕ�Ӧ�ľ���
Idx=knnsearch(pos,pos,'K',KNN);  %�ҵ����ڵĵ�idx������¼��¼
Idx = knnsearch(X,Y)  %��y����x����ڵ㣬һ������idx��
%% �ֶ�����a��c�ǶȺ�l1-l4
ang_a=[];
ang_c=[];
dis_l1=[];
dis_l2=[];
dis_l3=[];
dis_l4=[];
%%  ���ı������1743
S1=0.5.*sind(ang_a).*(0.024.*dis_l1).*(0.024.*dis_l4)  %��λ��nm2,,��������Ӧ����õ��(.*)
S2=0.5.*sind(ang_c).*(0.024.*dis_l2).*(0.024.*dis_l3)   %pxiel=0.024*0.024
S=S1+S2 
%% ���ı������1754
S1=0.5.*sind(ang_a).*(0.017.*dis_l1).*(0.017.*dis_l4)  %��λ��nm2,,��������Ӧ����õ��(.*)
S2=0.5.*sind(ang_c).*(0.017.*dis_l2).*(0.017.*dis_l3)   %pxiel=0.017*0.017
S=S1+S2 

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
M=27*1.6606e-27  % Al��ԭ������=���ԭ������*1.6606e-27 kg
E=80000    %��λkeV
m=9.1e-31   %��ֹ�ĵ���������������λkg
c=3e8    %����m/s
E1=2*M*E*(E+2*m*(c^2))
E2=((M+m)^2)*(c^2)+2*M*E
Em=E1/E2

%% 
M=72.6*1.6606e-27 %ԭ������
E=80000*1.602e-19 %J
m=9.1e-31   %��ֹ�ĵ���������������λkg
c=300000000    %����m/s

Em=((2*E*(E+2*m*(c^2)))/(M*(c^2)))/1.602e-19



