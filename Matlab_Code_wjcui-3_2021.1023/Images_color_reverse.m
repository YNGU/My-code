%% color reverse
I=imread('E:\++++++Matlab-Data\Exp data\Cu；Sn=1;1\ps\FFT HAADF 1448 20191219.tif'); %%imread(路径)阅读图片
figure,imshow(I); %%figure，imshow(变量)看图片
I_reverse=imcomplement(I); %%imcomplement(变量)反转图片
figure,imshow(I_reverse);
axis image
imwrite(rgbImage,'E:\++++++Matlab-Data\Exp data\Cu；Sn=1;1\ps\FFT HAADF 1613 20191219-gray.tif');
%% color reverse 2
I=imread('E:\Matlab-\cui\2019-1027-practice\16.14.15 CCD Acquire.tif');
figure,imshow(I);
%I_gray=rgb2gray(I); %%rgb2gray(变量）把rgb三通道图片变成灰度单通道图像
%figure,imshow(I_gray);
%I_reverse=imcomplement(I);
%figure,imshow(I_reverse);
I_reverse2=255-I; %%用灰度最高值255减去图片的灰度值达到反转的效果
figure,imshow(I_reverse2);

%% save mat_file
mat_file_name='E:\Matlab-\cui\2019-1027-practice\16.14.15 CCD Acquire.mat' %%保存mat文件，要加mat文件后缀
save(mat_file_name)
%% The axis transfomation of crystal orientation (3to4 axis)输入三轴参数，得到四轴晶向
clear
U=0;
V=-6;
W=1;
axis_3=[U,V,W];
u_1=(2*U-V)/3;
v_1=(2*V-U)/3;
t_1=-(u_1+v_1);
w_1=W;
axis_4_1=[u_1,v_1,t_1,w_1];%得到的形式为小数，应该继续弄个最小公倍数，转化为整数
%sym为符号变量存储方式，可以直接看到分数形式，之前为小数点/
%u_2=sym(u_1);
%v_2=sym(v_1);
%t_2=sym(t_1);
%w_2=sym(w_1);
gbs_1=lcm(u_1,v_1)
%% Q1:最小公倍数的计算公式
%% The axis transfamation (4to3 axis)四轴转变为三轴
clear
u=0;
v=1;
t=-1;
w=5;
axis_4=[u,v,t,w];
U=2*u+v;
V=2*v+u;
W=w;
axis_3=[U,V,W]
%% Orientation marix relationships（G_marix为E单位矩阵时为标准坐标轴rs）cos+acos函数结果为弧度制，cosd+acosd为角度制
clear %晶体坐标系与样品系统之间的关系用一个方向矩阵G来描述,在晶体系统中测量到的方向rc与样品系统测量到的方向rs的关系:rc=G*rs
alpha_1=acosd(0.707107);
beta_1=acosd(-0.408248);
gamma_1=acosd(0.577350);
alpha_2=acosd(0);
beta_2=acosd(0.816497);
gamma_2=acosd(0.577350);
alpha_3=acosd(-0.707107);
beta_3=acosd(-0.408248);
gamma_3=acosd(0.577350);
G_marix=[alpha_1,alpha_2,alpha_3;beta_1,beta_2,beta_3;gamma_1,gamma_2,gamma_3];
rs=[1 0 0; 0 1 0; 0 0 1];
rc=G_marix*rs %需要找到带轴与marix的关系，直接算出带轴方向
G_marix_correct=cosd(rc);%验证G_marix矩阵是否与VESTA里面的orientation marix相同
%% Q2:计算矩阵和带轴方向之间的关系 
a_deg=[alpha_1;alpha_2;alpha_3];
b_deg=[beta_1;beta_2;beta_3];
c_deg=[gamma_1;gamma_2;gamma_3];
UVW=cosd(a_deg)+cosd(b_deg)+cosd(c_deg);
%% 最小公倍数
