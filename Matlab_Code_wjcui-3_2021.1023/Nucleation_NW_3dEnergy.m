%% 先在matlab中取值，再将Z矩阵复制粘贴到origin中画图
clc
clear
%% 输入x，y数值，纳米液滴的形核公式
image = figure;  
d=0:0.00000000005:0.00000001;  %d为纳米线宽度
l=0:0.00000000005:0.00000001;  %l为纳米线长度
V=1; 
As=0.000000002; 
Ai1=0.00000000032; %以上三个值应该同ND形核时的参数设置一致
Ai2=0.00000000032; %Ai为γi2的设定值
[X,Y]=meshgrid(d,l); 
%for i=1:201;    %for循环去掉多余y>2x的值，前面x和y范围变了的话，这里也要变化
 %   for j=1:201;
  %      if Y(i,j)>2.*X(i,j);  %(i,j)表示矩阵中的对应第i行，第j列的元素
   %        Y(i,j)=NaN;  %NaN表示缺失值，没有数值
    %    end
    %end
%end
Gv=0.785.*V.*(X.^2).*Y; 
GAs=(3.14*As*X).*Y;         
GAi=0.785*(Ai1+Ai2).*(X.^2);  %两个界面的贡献
Z=-Gv+GAs+GAi;
surf(X,Y,Z);
shading interp  %插值使得图像变得平滑
title("ND ΔG energy")
xlabel('Width d'),ylabel('Length l'),zlabel('ΔG/ΔGv')
colormap(jet)
colorbar
hold on
savefig(image,'D:\++++++Data\Matlab Data\Exp_Data\NW_nucleation_jet')


