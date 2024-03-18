% in situ intensity analysis
clear
clc
%% 1. read the dm3 file  
dm3_file_name='D:\++++++Data\Matlab Data\Exp_Data\Test\1019 4.0 Mx HAADF STEM.dm3';

%% 2. read the dm3 file 
ImageSum_R=ReadDMFile(dm3_file_name);
ImageSum_R=double(ImageSum_R);

%% 3. display the image չ��ͼ��
imagesc(ImageSum_R);
axis image
axis off
colormap(inferno);
%caxis([0.02*256 0.3*256]);

%% 4. crop the image by defining the pixel area ��ͼ��pixel��ͼ
ImageSum_R_crop=ImageSum_R(:,470:595);   %�޸Ĳü��Ĳ���
h = fspecial('gaussian', [20 20],1);
ImageSum_BL=conv2(ImageSum_R_crop,h,'same');
ImageSum_BL(ImageSum_R_crop<0.0001)=0;
imagesc(ImageSum_BL);
axis image
axis off
colormap(inferno);
%% 4. crop the image if necessary  �ֶ�ѡ��ü������Ǳ���Ŀǰ����
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
%% 5. find local maxima % �ҵ�һ����Χ�ڵ��������ֵ,���ҵ�
BW = imregionalmax(ImageSum_BL); % imregionalmax�������ҳ���ǰ������Ԫ�ش��ڻ����ǿ��t��Ԫ�أ�������Ϊ1,�������ض���Ϊ0�������ض�ֵͼ��
[row,col] = find(BW); % �ҵ�BW�е�1�����ط�0��������row,colΪBW�����Ӧ���к���
%% 6. display the position of the peaks that wer found using maximum criteria %����ҵ���λ��
imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);
hold on  % ������ͼ�ϻ�
plot(col,row,'or','Markersize',4); % ��col��row�ĵ㣬plot(x,y)����x,y�������������Ǳ��������ͬ�ĳ��ȡ���������xΪ���ᣬ����y��

%% 7. use Gaussian distribution to refine atom column positions ͨ����˹��Ͼ�ϸ�ҵ㣬�ҵ�����7��������I��I0��x��y���ȣ�-��-
edge=10;  % ȥ����Ե��λ��
range=4;  % �Ե�ΪԲ�ģ��Ұ뾶Ϊ4*piexl��Բ�ڵĵ�����ϣ�
[fitresult_N,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,pos_refined]=...
    refine_atomic_columns(ImageSum_BL,[row col],edge,range,0);  % ������Ϊ0����ȫ��ͼ�ĵ㣬����100��100���㣿��һ��Ч��

%% 8. error check by compare image1 and image2  % ԭʼͼ�������ͼ�����������Ƿ���ȷ��ÿ�����λ��û�����Ե��������޹��򼴿�
figure
imagesc(image1-image2);  % Figure1��ԭʼ���ݣ�2���������
axis image

%% 9. calculate intensity of each peak  % 
row_gaussian=pos_refined(:,2);
col_gaussian=pos_refined(:,1);
pos_intensity=zeros(length(pos_refined),1); % ԭ����ǿ�����ݣ�pos_refinedΪ�Ż����ԭ����λ��
index=1;
for i=1:1:length(fitresult_N)
    if isempty(fitresult_N{i})
        continue;
    end
    pos_intensity(index)=fitresult_N{i}(1)+fitresult_N{i}(7);
    index=index+1;
end
%% 10. show if Gaussian distribution can improve the results ���һ�¾�ϸ��λ�öԲ���
figure
imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);
hold on
%plot(col,row,'or','Markersize',4);
plot(pos_refined(:,1),pos_refined(:,2),'ro','Markersize',4); 

%% 11. find 10 nearest neighbor of each point, and show an example ����10�����ڵ��λ�������
[Idx,D] = knnsearch(pos_refined,pos_refined,'K',10);
% calculate the angle between the peaks ͨ������Ƕ����ҽ��ڵ�
Angle=D;   % ���ڵ����ĽǶ�����
for i=1:1:length(pos_refined)
    x0=pos_refined(i,1);
    y0=pos_refined(i,2);
    for j=1:1:10
        x1=pos_refined(Idx(i,j),1);
        y1=pos_refined(Idx(i,j),2);
        Angle(i,j)=atan2(y1-y0,x1-x0)*180/pi();
    end
end

index=200;  % ������ֵλ�ö�Ӧ�ĵ�ĽǶ���������������øò������뿴�ĵ㣬����ͼfigure����ѡ�����Ĳ�����ʲô

imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);

hold on
plot(pos_refined(index,1),pos_refined(index,2),'xb','Markersize',4);
hold on
for i=1:1:10
    plot(pos_refined(Idx(index,i),1),pos_refined(Idx(index,i),2),'or','Markersize',4);
    text(pos_refined(Idx(index,i),1),pos_refined(Idx(index,i),2),num2str(Angle(index,i)));%�򿪽Ƕ���ͼ
    %text(pos_refined(Idx(index,i),1),pos_refined(Idx(index,i),2),num2str(D(index,i))); % �򿪾�����ͼ
    hold on
end

%% 12. use angle and distance pairs to find all those atoms ͨ���ǶȾ���ѡ��Ҫ�ĵ�
NNN=6;   % �Ĳ�������Ҫ�������ڵ�
NN_angles=[-1 55 145 180 -124 -35];   % ��Ҫ��ȡ���ڵ�Ĵ��½Ƕ�
NN_dist=[16 9.5 13 16.5 8.9 13.5];   % ��Ҫ��ȡ���ڵ�Ĵ�֮���룬����һ����%��������ʱ��˳��д����
angle_threshold=20;    % �Ƕȵ���Χ������
dist_threshold=5;    % �������Χ��pixel
pos_NN=zeros(length(pos_refined),NNN);    % ���ڵ������
pos_angle=zeros(length(pos_refined),NNN);   % ��ȡ���ڵ�ĽǶ�
pos_dis=zeros(length(pos_refined),NNN);   % ��ȡ���Ľ��ڵ�ľ���
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

%% 13. plot peaks that can find six neareast neighbors ����ҵĵ�Բ���
figure
imagesc(ImageSum_BL);
axis image
axis off
colormap(gray);
hold on
%plot(col,row,'or','Markersize',4);
x1=col_gaussian(pos_NNN==4);
y1=row_gaussian(pos_NNN==4);
plot(x1',y1','bo','Markersize',4);  % ��ɫ�����Χ���ĸ���
x1=col_gaussian(pos_NNN==5);
y1=row_gaussian(pos_NNN==5);
plot(x1',y1','ro','Markersize',4);  % ��ɫ��������
x1=col_gaussian(pos_NNN==6);
y1=row_gaussian(pos_NNN==6);
plot(x1',y1','yo','Markersize',4);  % ��ɫ������������ڵ�ĵ�

%% 14. index all the peaks 
initial_peak=1200;  % �Ĳ������Բ�������һ���ĵ㶼����
peak_index=nan(length(pos_refined),2);
peak_index(initial_peak,1)=1;
peak_index(initial_peak,2)=1;
peak_index=index_single_peak_1(pos_NN,pos_NNN,initial_peak, peak_index,5,1,2,4);  % ��ǽ��ڵ�ķ�����index
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
%% 15. show the index of all the peaks  ͼ�п�index=[i, j]�Բ���
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
%% 16. draw distance distribution along [1 0] ����[1 0] ���򻭾���
[dist_matrix1,f]=draw_distance_between_peaks(index_matrix(:,:),pos_refined,ImageSum_BL,1,6,16,8,13,2,4,[1 0]); %�Ĳ���pixel�����������Сֵ*2����ɫ��Χ*2���߿�Բ�δ�С��[����]
%xlim([250 450]);  % dist_matrixl ����������ݣ�draw_distance_between_peaks����������
%ylim([180 380]);
%print(f,'-dtiff', '-r600', [dm3_file_name(1:end-4),'_figure1.tiff']);
%print(f,'-dpng',  '-r600', [dm3_file_name(1:end-4),'figure1.png']); % ����50*50cmͼƬ,��draw_distance���������

%% 17. draw distance distribution along [0 1]
[dist_matrix2,f]=draw_distance_between_peaks(index_matrix(:,:),pos_refined,ImageSum_BL,1,15,17.5,0,0,2,4,[0 1]);
%xlim([250 450]);
%ylim([180 380]);
%print(f,'-dtiff', '-r600', [dm3_file_name(1:end-4),'_figure2.tiff']);
%print(f,'-dpng',  '-r600', [dm3_file_name(1:end-4),'figure2.png']);
%% 18. find reasonalbe upper and lower bounds for distance plot ����ֱ��ͼ
hist(dist_matrix2(:),0:0.1:20);   %ͳ��dist_matrix2�ľ������ݣ����Ը�dist_matrix
%% 19. find reasonable upper and lower bounds for intensity plot ��intensityֱ��ͼ
hist(pos_intensity(:));  
%% 20. draw rectangles showing the intensity of the peaks ����intensit mapͼ
[int_matrix,f]=intensity_mapping(index_matrix(:,:),pos_refined,ImageSum_BL,pos_intensity,0,0,150,250,0.7,7,2); % ǿ�ȷ�Χ���߿����
%xlim([250 450]);
%ylim([180 380]);
%print(f,'-dpng',  '-r600', [dm3_file_name(1:end-4),'figure3.png']); 
%% 21. save the mat file
mat_file_name=[dm3_file_name(1:end-3) 'mat'];
save(mat_file_name);