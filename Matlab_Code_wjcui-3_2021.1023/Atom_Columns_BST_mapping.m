% in situ intensity analysis
clear
clc
%% 1. read the dm3 file  
dm3_file_name='D:\++++++Data\Matlab Data\Exp_Data\Test\1019 4.0 Mx HAADF STEM.dm3';

%% 2. read the dm3 file 
ImageSum_R=ReadDMFile(dm3_file_name);
ImageSum_R=double(ImageSum_R);

%% 3. display the image 展现图像
imagesc(ImageSum_R);
axis image
axis off
colormap(inferno);
%caxis([0.02*256 0.3*256]);

%% 4. crop the image by defining the pixel area 按图像pixel截图
ImageSum_R_crop=ImageSum_R(:,470:595);   %修改裁剪的参数
h = fspecial('gaussian', [20 20],1);
ImageSum_BL=conv2(ImageSum_R_crop,h,'same');
ImageSum_BL(ImageSum_R_crop<0.0001)=0;
imagesc(ImageSum_BL);
axis image
axis off
colormap(inferno);
%% 4. crop the image if necessary  手动选择裁剪，但是变量目前不对
figure;
imagesc(ImageSum_R);
axis image;
colormap(inferno);

h = imrect;

position = wait(h);
position=floor(position);
image_rangex=position(2):1:position(2)+position(4);
image_rangey=position(1):1:position(1)+position(3);
ImageSum_R=ImageSum_R(image_rangex,image_rangey);
imagesc(ImageSum_R);
axis image
%% 5. find local maxima % 找到一定范围内的像素最大值,粗找点
BW = imregionalmax(ImageSum_BL); % imregionalmax函数：找出当前区域内元素大于或等于强度t的元素，将其置为1,其他像素都设为0。并返回二值图像。
[row,col] = find(BW); % 找到BW中的1，返回非0的索引，row,col为BW里面对应的行和列
%% 6. display the position of the peaks that wer found using maximum criteria %标记找到的位置
imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);
hold on  % 保持在图上画
plot(col,row,'or','Markersize',4); % 画col，row的点，plot(x,y)中若x,y是向量，则它们必须具有相同的长度。函数将以x为横轴，绘制y。

%% 7. use Gaussian distribution to refine atom column positions 通过高斯拟合精细找点，找点至少7个参数，I，I0，x，y，θ，-，-
edge=10;  % 去掉边缘的位置
range=4;  % 以点为圆心，找半径为4*piexl的圆内的点做拟合，
[fitresult_N,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,pos_refined]=...
    refine_atomic_columns(ImageSum_BL,[row col],edge,range,0);  % 改数字为0是找全幅图的点，先用100找100个点？试一试效果

%% 8. error check by compare image1 and image2  % 原始图像与拟合图像相减看拟合是否正确，每个拟合位置没有明显的线条，无规则即可
figure
imagesc(image1-image2);  % Figure1是原始数据，2是拟合数据
axis image

%% 9. calculate intensity of each peak  % 
row_gaussian=pos_refined(:,2);
col_gaussian=pos_refined(:,1);
pos_intensity=zeros(length(pos_refined),1); % 原子列强度数据，pos_refined为优化后的原子列位置
index=1;
for i=1:1:length(fitresult_N)
    if isempty(fitresult_N{i})
        continue;
    end
    pos_intensity(index)=fitresult_N{i}(1)+fitresult_N{i}(7);
    index=index+1;
end
%% 10. show if Gaussian distribution can improve the results 检查一下精细点位置对不对
figure
imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);
hold on
%plot(col,row,'or','Markersize',4);
plot(pos_refined(:,1),pos_refined(:,2),'ro','Markersize',4); 

%% 11. find 10 nearest neighbor of each point, and show an example 计算10个近邻点的位置与距离
[Idx,D] = knnsearch(pos_refined,pos_refined,'K',10);
% calculate the angle between the peaks 通过计算角度来找近邻点
Angle=D;   % 近邻点计算的角度数据
for i=1:1:length(pos_refined)
    x0=pos_refined(i,1);
    y0=pos_refined(i,2);
    for j=1:1:10
        x1=pos_refined(Idx(i,j),1);
        y1=pos_refined(Idx(i,j),2);
        Angle(i,j)=atan2(y1-y0,x1-x0)*180/pi();
    end
end

index=200;  % 看该数值位置对应的点的角度情况（度数），该该参数看想看的点，并画图figure，来选择后面的参数用什么

imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);

hold on
plot(pos_refined(index,1),pos_refined(index,2),'xb','Markersize',4);
hold on
for i=1:1:10
    plot(pos_refined(Idx(index,i),1),pos_refined(Idx(index,i),2),'or','Markersize',4);
    text(pos_refined(Idx(index,i),1),pos_refined(Idx(index,i),2),num2str(Angle(index,i)));%打开角度视图
    %text(pos_refined(Idx(index,i),1),pos_refined(Idx(index,i),2),num2str(D(index,i))); % 打开距离视图
    hold on
end

%% 12. use angle and distance pairs to find all those atoms 通过角度距离选需要的点
NNN=6;   % 改参数，需要几个近邻点
NN_angles=[-1 55 145 180 -124 -35];   % 需要获取近邻点的大致角度
NN_dist=[16 9.5 13 16.5 8.9 13.5];   % 需要获取近邻点的大之距离，从上一步的%来看，逆时针顺序写参数
angle_threshold=20;    % 角度的误差范围，度数
dist_threshold=5;    % 距离的误差范围，pixel
pos_NN=zeros(length(pos_refined),NNN);    % 近邻点的数量
pos_angle=zeros(length(pos_refined),NNN);   % 获取近邻点的角度
pos_dis=zeros(length(pos_refined),NNN);   % 获取到的近邻点的距离
pos_NNN=zeros(length(pos_refined),1);   
for i=1:1:length(pos_refined)
    for j=1:1:NNN
        
        for k=1:1:10
            angle_diff=abs(Angle(i,k)-NN_angles(j));
            
            if abs(D(i,k)-NN_dist(j))<=dist_threshold
                if angle_diff<=angle_threshold || abs(angle_diff-360)<=angle_threshold
                    pos_NN(i,j)=Idx(i,k);
                    pos_angle(i,j)=Angle(i,k);
                    pos_dis(i,j)=D(i,k);
                    pos_NNN(i)=pos_NNN(i)+1;
                end
            end
        end
    end
end

%% 13. plot peaks that can find six neareast neighbors 检查找的点对不对
figure
imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);
hold on
%plot(col,row,'or','Markersize',4);
x1=col_gaussian(pos_NNN==4);
y1=row_gaussian(pos_NNN==4);
plot(x1',y1','bo','Markersize',4);  % 蓝色标记周围有四个点
x1=col_gaussian(pos_NNN==5);
y1=row_gaussian(pos_NNN==5);
plot(x1',y1','ro','Markersize',4);  % 红色标记五个点
x1=col_gaussian(pos_NNN==6);
y1=row_gaussian(pos_NNN==6);
plot(x1',y1','yo','Markersize',4);  % 黄色标记有六个近邻点的点

%% 14. index all the peaks 
initial_peak=1200;  % 改参数，试参数让下一步的点都覆盖
peak_index=nan(length(pos_refined),2);
peak_index(initial_peak,1)=1;
peak_index(initial_peak,2)=1;
peak_index=index_single_peak_1(pos_NN,pos_NNN,initial_peak, peak_index,5,1,2,4);  % 标记近邻点的方向与index
% 5 -> [1 0]
% 1 -> [0 1]
% 2 -> [-1 0]
% 4 -> [0 -1]
peak_index(:,1)=peak_index(:,1)-min(peak_index(:,1))+1;
peak_index(:,2)=peak_index(:,2)-min(peak_index(:,2))+1;
ps=max(peak_index);
index_matrix=zeros(ps);
for i=1:1:length(pos_refined)
    if ~isnan(peak_index(i,1))
        index_matrix(peak_index(i,1),peak_index(i,2))=i;
    end
end
imagesc(index_matrix);
%% 15. show the index of all the peaks  图中看index=[i, j]对不对
imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);

hold on

for i=1:1:length(pos_refined)
    if isnan(peak_index(i,1))   
        continue;
    end
    plot(pos_refined(i,1),pos_refined(i,2),'or','Markersize',4);
    text(pos_refined(i,1),pos_refined(i,2),...
        [num2str(peak_index(i,1)) ',' num2str(peak_index(i,2))]);
    hold on
end
%% 16. draw distance distribution along [1 0] 沿着[1 0] 方向画距离
[dist_matrix1,f]=draw_distance_between_peaks(index_matrix(:,:),pos_refined,ImageSum_BL,1,6,16,8,13,2,4,[1 0]); %改参数pixel，距离最大最小值*2，颜色范围*2，线宽，圆形大小，[方向]
%xlim([250 450]);  % dist_matrixl 储存距离数据，draw_distance_between_peaks看参数设置
%ylim([180 380]);
%print(f,'-dtiff', '-r600', [dm3_file_name(1:end-4),'_figure1.tiff']);
%print(f,'-dpng',  '-r600', [dm3_file_name(1:end-4),'figure1.png']); % 保存50*50cm图片,在draw_distance函数里面改

%% 17. draw distance distribution along [0 1]
[dist_matrix2,f]=draw_distance_between_peaks(index_matrix(:,:),pos_refined,ImageSum_BL,1,15,17.5,0,0,2,4,[0 1]);
%xlim([250 450]);
%ylim([180 380]);
%print(f,'-dtiff', '-r600', [dm3_file_name(1:end-4),'_figure2.tiff']);
%print(f,'-dpng',  '-r600', [dm3_file_name(1:end-4),'figure2.png']);
%% 18. find reasonalbe upper and lower bounds for distance plot 距离直方图
hist(dist_matrix2(:),0:0.1:20);   %统计dist_matrix2的距离数据，可以改dist_matrix
%% 19. find reasonable upper and lower bounds for intensity plot 画intensity直方图
hist(pos_intensity(:));  
%% 20. draw rectangles showing the intensity of the peaks 导出intensit map图
[int_matrix,f]=intensity_mapping(index_matrix(:,:),pos_refined,ImageSum_BL,pos_intensity,0,0,150,250,0.7,7,2); % 强度范围，线宽参数
%xlim([250 450]);
%ylim([180 380]);
%print(f,'-dpng',  '-r600', [dm3_file_name(1:end-4),'figure3.png']); 
%% 21. save the mat file
mat_file_name=[dm3_file_name(1:end-3) 'mat'];
save(mat_file_name);