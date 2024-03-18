%% color reverse
I=imread('E:\++++++Matlab-Data\Exp data\Cu��Sn=1;1\ps\FFT HAADF 1448 20191219.tif'); %%imread(·��)�Ķ�ͼƬ
figure,imshow(I); %%figure��imshow(����)��ͼƬ
I_reverse=imcomplement(I); %%imcomplement(����)��תͼƬ
figure,imshow(I_reverse);
axis image
imwrite(rgbImage,'E:\++++++Matlab-Data\Exp data\Cu��Sn=1;1\ps\FFT HAADF 1613 20191219-gray.tif');
%% color reverse 2
I=imread('E:\Matlab-\cui\2019-1027-practice\16.14.15 CCD Acquire.tif');
figure,imshow(I);
%I_gray=rgb2gray(I); %%rgb2gray(��������rgb��ͨ��ͼƬ��ɻҶȵ�ͨ��ͼ��
%figure,imshow(I_gray);
%I_reverse=imcomplement(I);
%figure,imshow(I_reverse);
I_reverse2=255-I; %%�ûҶ����ֵ255��ȥͼƬ�ĻҶ�ֵ�ﵽ��ת��Ч��
figure,imshow(I_reverse2);

%% save mat_file
mat_file_name='E:\Matlab-\cui\2019-1027-practice\16.14.15 CCD Acquire.mat' %%����mat�ļ���Ҫ��mat�ļ���׺
save(mat_file_name)
%% The axis transfomation of crystal orientation (3to4 axis)��������������õ����ᾧ��
clear
U=0;
V=-6;
W=1;
axis_3=[U,V,W];
u_1=(2*U-V)/3;
v_1=(2*V-U)/3;
t_1=-(u_1+v_1);
w_1=W;
axis_4_1=[u_1,v_1,t_1,w_1];%�õ�����ʽΪС����Ӧ�ü���Ū����С��������ת��Ϊ����
%symΪ���ű����洢��ʽ������ֱ�ӿ���������ʽ��֮ǰΪС����/
%u_2=sym(u_1);
%v_2=sym(v_1);
%t_2=sym(t_1);
%w_2=sym(w_1);
gbs_1=lcm(u_1,v_1)
%% Q1:��С�������ļ��㹫ʽ
%% The axis transfamation (4to3 axis)����ת��Ϊ����
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
%% Orientation marix relationships��G_marixΪE��λ����ʱΪ��׼������rs��cos+acos�������Ϊ�����ƣ�cosd+acosdΪ�Ƕ���
clear %��������ϵ����Ʒϵͳ֮��Ĺ�ϵ��һ���������G������,�ھ���ϵͳ�в������ķ���rc����Ʒϵͳ�������ķ���rs�Ĺ�ϵ:rc=G*rs
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
rc=G_marix*rs %��Ҫ�ҵ�������marix�Ĺ�ϵ��ֱ��������᷽��
G_marix_correct=cosd(rc);%��֤G_marix�����Ƿ���VESTA�����orientation marix��ͬ
%% Q2:�������ʹ��᷽��֮��Ĺ�ϵ 
a_deg=[alpha_1;alpha_2;alpha_3];
b_deg=[beta_1;beta_2;beta_3];
c_deg=[gamma_1;gamma_2;gamma_3];
UVW=cosd(a_deg)+cosd(b_deg)+cosd(c_deg);
%% ��С������
