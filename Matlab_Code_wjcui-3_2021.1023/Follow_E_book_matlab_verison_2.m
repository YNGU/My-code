%% 2021.1211  ͼ������������jpg��tif�ĸ�ʽ����
image=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH01_Images\Fig0104.tif.');
imshow(image);
%figure, imshow(g);    %��figure������һ��ͼ������ʾ�ڶ���ͼ
%imwrite(image,'D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH01_Images\1.jpg','quality',50) ...
    %��ͼ��λjpg��ʽ�����ϡ�quality�������������ֵ0-100��ʾͼ��������ԽСͼ��Խ�ֲ�
imwrite(image,'D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH01_Images\1.tif', 'compression', ...
    'none', 'resolution', [1024 1024])     %����tif�ļ��ı�׼��ʽ��'compression', 'parameter', 'resolution', [rolres rowres]�� parameter��com��ѹ��������colres��rowres�������еķֱ���
%% 2021.1213   ��������
v=[1 3 5 7 9];
v_2=v(2);   %()��ʾ������λ��,��������С���ű�ʾ����
w=v.';    %ת�þ���
v(1:3);    %:��ʹ��
v(1:end);  %end�������һ����ֵ
v([1 3 5]);  %��Ӧ��������
v(1:2:end);    %���ò���
A=[1 2 3; 4 5 6; 7 8 9];
%A(2,:);   %:��ʾһ���е�����Ԫ��
A(2,1:end)  %����һ�е�Ч�� 
D=logical([1 0 0; 0 0 1; 0 0 0]);
A(D);
%% �������
f=@sin;   %@�������ú������ --��һ�ຯ�����-�򵥺������
a=f(pi/4);
g=@(x) x.^2     % �ڶ����������-�����������==��һ����ʽ��@��input argument list�� expression
g(2);
r=@(x,y) sqrt(x.^2 + y.^2);
r(3, 4);
%% 2021.1214 imadjust�Ҷ�ת��
image=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH02_Images\Fig0203(a).tif');
image_contrast1=imadjust(image,[0 1],[1 0],1);    %imadjust������out�������Сֵ�෴���ɻҶȷ�ת��ֵ������0-1֮��
image_contrast2=imcomplement(image);              %imcomplementҲ�ǻҶȷ�ת
figure(2);
subplot(2,2,1),imshow(image);
subplot(2,2,2),imshow(image_contrast1);
subplot(2,2,3),imshow(image_contrast2);
image2=imadjust(image,[0.5 0.75],[0 1],1);
subplot(2,2,4),imshow(image2);
image3=imadjust(image,[0 1],[0.25 0.75],1);
figure,imshow(image3);                         %Ҫ�����һ��ͼƬ���ڣ�ǰ���figure��
%% ����gammaֵ����1С��1ʱ�ı仯��Ĭ��ֵΪ1
image=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH02_Images\Fig0203(a).tif')
image_tran1=imadjust(image,[0 1],[0.2 0.9],1)
image_tran2=imadjust(image,[0 1],[0.2 0.9],2)
image_tran3=imadjust(image,[0 1],[0.2 0.9],3)
image_tran4=imadjust(image,[0 1],[0.2 0.9],0.7)
image_tran5=imadjust(image,[0 1],[0.2 0.9],0.3)
image_tran6=imadjust(image,[0 1],[0.2 0.9],0)      %gammaΪ0ʱ��û��ͼ����
figure(1)
subplot(3,2,1),imshow(image_tran1);
subplot(3,2,2),imshow(image_tran2);
subplot(3,2,3),imshow(image_tran3);
subplot(3,2,4),imshow(image_tran4);
subplot(3,2,5),imshow(image_tran5);
subplot(3,2,6),imshow(image_tran6);
%% stretchlim(f)����һ����Ԫ����������һ�����޺͸������
%% log�Ա���Աȶ�����任������һ����ֵ���ڸ��˵�ֵ���źţ�Ϊ�˽���ֵ���ź�Ҳ���ֳ���
f=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH02_Images\Fig0205(a).tif')
%imshow(f)
f1=mat2gray(f);                         %mat2gray��ͼƬ��pixelֵ�޶���[0,1]֮��
f2=im2uint8(f1);                           %im2uint8��ͼƬ��Ϊunit8��ʽ���ң�pixelֵ�޶���[0,255]֮��
g=im2uint8(mat2gray(log(1+double(f))));      %log=loge��log2=log2
figure;
subplot(2,1,1),imshow(f);
subplot(2,1,2),imshow(g);
%% С���ɳ�Ϊ����������Ϊһ��С��3.14�����ʹ��ָ��������ʽ�Ļ���3.14E0��С����������⸡��
%% 2021.1216 ͼ��ֱ��ͼ+++���Ϸֺ�֮��Ͳ��������������д������У�ֱ�Ӽ���ʡ��ʱ�䣬
f=imread('D:\++++++Data\Matlab Data\e_Book\DIPUM2E_International_Original_Book_Images\DIPUM2E_International_Ed_CH02_Images\Fig0205(a).tif'); 
h=imhist(f,25);     %imhist(f,N) N�ǰ�����������Ϊ���㣬���������N
horz=linspace(0,255,25);    %linspace(x1,x2,N)��x1��x2��N�ֱ�Ϊ��ʼֵ����ֵֹ��Ԫ�صĸ�����NĬ�ϵ���Ϊ100��
bar(horz,h);       %bar,stem,plot�����������Ļ���
stem(horz,h);
plot(horz,h);
%imhist(f);   %imhistֱ�ӿ���չʾͼƬ�ĻҶȷֲ�
%figure(1);
%subplot(2,2,1),imhist(f);
%subplot(2,2,2),bar(horz,h);
%subplot(2,2,3),stem(horz,h);
%subplot(2,2,4),plot(horz,h);     %����ֱ����subplot������ź���������imshow����Ҳ���ԡ�
%%  xy��ķ�Χ����������
imhist(f);
axis([0 100 0 60000]);               %axis����xy�����Сֵ�����ֵ
set(gca, 'xtick', 0:10:100);
set(gca, 'ytick', 0:20000:60000);      % ����x��y��ļ����С
xlabel('x axis', 'fontsize', 20);
ylabel('y axis', 'fontsize', 20);
title('gray distubrition','fontsize',20);   %������Ŀ
ylim('auto');
xlim('auto');

%%

