%% 先在matlab中取值，再将Z矩阵复制粘贴到origin中画图
clc
clear
%% 输入x，y数值，纳米液滴的形核公式
image = figure;  %新建一个图像
%set(gcf,'unit','centimeters','position',[0 0 6 4]); %设置画布大小，单位厘米
%set(gcf,'color',[1  1]) %设置背景颜色，x/255=？，范围为0-1
r=0:0.00000000005:0.000000008;  %单位m
h=0:0.00000000005:0.000000008;  
V=1; %V为ΔGv/ΔGv
As=0.000000002; %As为γs/ΔGv
Ai=0.00000000032; %Ai为γi/ΔGv
[X,Y]=meshgrid(r,h); 
for i=1:161;    %for循环去掉多余y>2x的值，前面x和y范围变了的话，这里也要变化
    for j=1:161;
        if Y(i,j)>2.*X(i,j);  %(i,j)表示矩阵中的对应第i行，第j列的元素
           Y(i,j)=NaN;  %NaN表示缺失值
        end
    end
end
Gv=1.047.*V.*(Y.^2).*[(3*X)-Y];  %体能贡献 对过origin的图，没问题
GAs=(6.28*As*X).*Y;   %表面能贡献
GAi=3.14*Ai.*(2.*Y.*X-Y.^2)   %(X.^2).*[1-((Y-X).^2)./(X.^2)];
%.^2就是矩阵中每个元素的平方,矩阵对应元素相乘用.*,./也一样
Z=-Gv;
surf(X,Y,Z);
shading interp  %插值使得图像变得平滑
title("ND ΔG/ΔGv energy")
xlabel('Radius r'),ylabel('Height h'),zlabel('ΔG/ΔGv')
colormap(jet)
colorbar

%zticks([])；%让z轴刻度不显示
hold on
savefig(image,'D:\++++++Data\Matlab Data\Exp_Data\ND_3dEnergy')


