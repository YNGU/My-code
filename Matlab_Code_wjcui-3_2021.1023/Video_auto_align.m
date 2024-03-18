%% clear all
clear
clc
%% 1. assign the video file name
video_file_name='D:\++++++Data\Matlab Data\Exp_Data\test\2042 4.0 Mx 20231002 STEM.mrc';
%% set mat file name
mat_file_name=[video_file_name(1:end-3) 'mat'];

%% 2. get mrc data 给video文件转换成double格式
data_mrc = ReadMRC(video_file_name);
data_mrc = double(data_mrc);
%% 3. the number of frames
s1=size(data_mrc);
nof=s1(3);

%% 4. show the first frame

figure;
imagesc(data_mrc(:,:,1));
axis image
colormap(gray);

%% 没有运行的步骤
figure;
imagesc(mov(1).cdata);
axis image;
colormap(hot);

h = imrect;

position = wait(h);
position=floor(position);
image_rangex=position(2):1:position(2)+position(4);
image_rangey=position(1):1:position(1)+position(3);
image_area=mov(1).cdata(image_rangex,image_rangey);
imagesc(image_area);
axis image
%% 5. 
[folder,name,ext] = fileparts(video_file_name);

%% 6. 把视频变成图片放在一个文件夹里面，key_frames文件夹
if ~exist([folder,'/key_frames'], 'dir')
    mkdir([folder,'/key_frames']); 
end

for k=1:1:nof
    frame_k=data_mrc(:,:,k);
    f_max=max(frame_k(:));
    f_min=min(frame_k(:)); %减去图像衬度最小数值
    imwrite((data_mrc(:,:,k)-f_min)/(f_max-f_min),[folder,'/key_frames/frame'...
        ,num2str(k),'.jpg'],'jpg');
    ImageFrames{k}=(data_mrc(:,:,k)-f_min)/(f_max-f_min)*255; %255改大，图像衬度整体变亮，暗的地方能看清
end

%% 7. 给图片加一个blur的滤波
% 1 for 16 nm X 16 nm
% 2 for 8 nm X 8 nm
blur_index=1;
%% 8. now save all the frames differently, 使图像变色，选择合适的参数
if ~exist([folder,'/key_frames_blurred'], 'dir')
    mkdir([folder,'/key_frames_blurred']);
end


h = fspecial('gaussian', [20 20],blur_index);

for i=1:1:length(ImageFrames)
    ImageSum_BL=conv2(ImageFrames{i},h,'same');
    J=imwrite_color(ImageSum_BL,10,255); %10,200为最大最小阈值，改变两者参数调节图片看哪个更美观
    imwrite(J,[folder,'/key_frames_blurred/frame'...
            ,num2str(i),'.jpg'],'jpg');
end

%% 可以试一下，看哪个好一些，和上一步同功能，算法不同，
if ~exist([folder,'\denoise'], 'dir')
    mkdir([folder,'\denoise']);
end
net = denoisingNetwork('DnCNN');


h = fspecial('gaussian', [20 20],2);

for i=1:1:length(ImageFrames)    
    denoisedI = denoiseImage(ImageFrames{i},net);
    ImageSum_BL=conv2(denoisedI,h,'same');
    
    J=imwrite_color(ImageSum_BL,10,200);%10,200为最大最小阈值，改变两者参数调节图片看哪个更美观
    
    imwrite( J,[folder,'\denoise\frame'...
        ,num2str(i),'.jpg'],'jpg');
end
%% 9. exclude some frames 排除一些无用的图片
exclude_frames=[99:1:147];
frames=1:1:length(ImageFrames);
frames(exclude_frames)=[];
% frames=1:1:57;
%% 10 align the images to its previous one 给图片align
if ~exist([folder,'/aligned'], 'dir')
    mkdir([folder,'/aligned']);
end

first_frame=ImageFrames{frames(1)};

f=figure;
imagesc(first_frame);
axis image;
colormap(gray);

h = imrect;
position = wait(h);
position=floor(position);
cc_rangex=position(2):1:position(2)+position(4);
cc_rangey=position(1):1:position(1)+position(3);
close(f);

shifts=zeros(length(frames),2);

for i=2:1:length(frames)
    first_frame=ImageFrames{frames(i-1)};
    current_frame=ImageFrames{frames(i)};
    
    data1=clear_background(first_frame,10,0);
    data2=clear_background(current_frame,10,0);
    data1=data1(cc_rangex,cc_rangey);
    data2=data2(cc_rangex,cc_rangey);
    
    s=dftregistration((fft2(data1)),(fft2(data2)),10);
    
    [i s(3) s(4)]
    shifts(i,1)=s(3)+shifts(i-1,1);
    shifts(i,2)=s(4)+shifts(i-1,2);
end
shifts_adjacent=zeros(length(frames),2);
for i=2:1:length(frames)
    shifts_adjacent(i,1)=shifts(i,1)-shifts(i-1,1);
    shifts_adjacent(i,2)=shifts(i,2)-shifts(i-1,2);
end
%% 11. save the images to compare if the shifts were calculated correctly 对比前后一张图片
if ~exist([folder,'/check_cc'], 'dir')
    mkdir([folder,'/check_cc']);
end

for k=2:1:length(frames)
    Image_shifted=imtranslate(ImageFrames{frames(k)},[shifts_adjacent(k,2) shifts_adjacent(k,1)]);
    ImageSum_BL=[ImageFrames{frames(k-1)} Image_shifted];
    J = imadjust(ImageSum_BL/255,[0;1], [0;1]);
    imwrite(J,[folder,'/check_cc/frame'...
            ,num2str(k),'.jpg'],'jpg');
end
%% this give an opportunity to fix the misalignment 可以不运行的步骤
% click on the first image and then the second image
for k=1:1:20
    subplot(1,3,1);
    imagesc(ImageFrames{frames(k-1)});
    axis image
    colormap(gray);
    title(['x=' num2str(shifts_adjacent(k,1)) ' y=' num2str(shifts_adjacent(k,2))]);
    subplot(1,3,2);
    imagesc(ImageFrames{frames(k)});
    axis image
    colormap(gray);
    subplot(1,3,3);
    imagesc(imtranslate(ImageFrames{frames(k)},[shifts_adjacent(k,2) shifts_adjacent(k,1)]));
    axis image
    colormap(gray);
    [x1, y1, but] = ginput(1);
    if but == 110 % 'N'
        continue;
    end
    [x2, y2, but] = ginput(1);
    title(['x=' num2str(y1-y2) ' y=' num2str(x1-x2)]);
    shifts_adjacent(k,1)=y1-y2;
    shifts_adjacent(k,2)=x1-x2;
    [x2, y2, but] = ginput(1);
end
%% 12. restore shifts 保存偏移量
for i=2:1:length(frames)
    shifts(i,1)=shifts_adjacent(i,1)+shifts(i-1,1);
    shifts(i,2)=shifts_adjacent(i,2)+shifts(i-1,2);
end
%% 13. save the aligned images 最后一步，保存到align文件夹，已经align好的图片导出
x_left =ceil(abs(min(shifts(:,1))))+1;
x_right=ceil(abs(max(shifts(:,1))))+1;
y_top   =ceil(abs(min(shifts(:,2))))+1;
y_bottom=ceil(abs(max(shifts(:,2))))+1;
[x_left x_right y_top y_bottom]
[sizex sizey]=size(ImageFrames{1});
image_till=sizex; % this deals with the time tag
for i=1:1:length(frames)
    large_frame=zeros(x_left+sizex+x_right,y_top+sizey+y_bottom);
    %large_frame(x_left:x_left+image_till-1,y_top:y_top+sizey-1)=...
    %    ImageFrames{frames(i)}(1:image_till,:);
    large_frame(x_left:x_left+image_till-1,y_top:y_top+sizey-1)=...
        ImageFrames{frames(i)}(1:image_till,:);
    Image_coarse_align{i}=imtranslate(large_frame,[shifts(i,2) shifts(i,1)]);
    %     current_frame(x_left:x_left+sizex-1,y_top+image_till:y_top+sizey)=...
    %         ImageFrames{i}(:,image_till:sizey);
    h = fspecial('gaussian', [20 20],1);
    %ImageSum_BL=conv2(Image_coarse_align{i},h,'same');
    ImageSum_BL=Image_coarse_align{i};
    %J = imadjust(Image_coarse_align{i}/255,[0;1], [0;1]);
    J=imwrite_color(ImageSum_BL,10,255);
    imwrite(J,[folder,'/aligned/frame'...
            ,num2str(i),'.jpg'],'jpg');
end
%% check the average of the frames
average_intensity=zeros(length(frames),1);
for i=1:1:length(frames)
    temp=ImageFrames{frames(i)}(1:image_till,:);
    average_intensity(i)=mean(temp(:));
    max_intensity(i)=max(temp(:));
    min_intensity(i)=min(temp(:));
end
figure
hold all
plot(average_intensity,'o');
plot(min_intensity);
plot(max_intensity);
%% crop it to view intesting areas
if ~exist([folder,'\roi'], 'dir')
    mkdir([folder,'\roi']);
end
first_frame=Image_coarse_align{29};
% for i=75:2:75
% first_frame=first_frame+... 
%     Image_coarse_align{i};
% end
f=figure;
imagesc(first_frame);
axis image;
colormap(gray);

h = imrect;
position = wait(h);
position=floor(position);
cc_rangex=position(2):1:position(2)+position(4);
cc_rangey=position(1):1:position(1)+position(3);
close(f);

for i=1:1:length(frames)
    current_frame=...
        Image_coarse_align{i}(cc_rangex,cc_rangey);
    %h = fspecial('gaussian', [20 20],blur_index);
    %ImageSum_BL=conv2(current_frame,h,'same');
    ImageSum_BL=current_frame;
    %J = imadjust(ImageSum_BL/255,[0;average_intensity(i)*3/255], [0;1]);
    J = ImageSum_BL;%imadjust(ImageSum_BL/255,[0;1], [1;255]);
    J=imwrite_color(J,1,255);
    time_tag=round((key_frame_list(i)-key_frame_list(1))/12);
%     roi_list{i} = insertText(J,[0 0],['t=' num2str(time_tag) 's'],'FontSize',48,'BoxColor',...
%     'black','BoxOpacity',1,'TextColor','white');
    roi_list{i}=J;
    imwrite( roi_list{i},[folder,'\roi\frame'...
            ,num2str(i),'.jpg'],'jpg');
end


%% remove the background
N=10;
ImageSum_f=ordfilt2(ImageSum_R_blurred,1,ones(N,N));
Image_index=(ImageSum_R_blurred-ImageSum_f)';
Image_index=2*Image_index/max(Image_index(:));
Image_index=255*Image_index;
Image_index(Image_index<0)=0;
Image_index(Image_index>255)=255;
Image_index=round(Image_index);
Image_index=Image_index+1;
ImageSum_RGB = ind2rgb(Image_index,inferno(256));
imshow(ImageSum_RGB);
png_file_name=[dm3_file_name(1:end-4) '_bkgd.png'];
imwrite(ImageSum_RGB,png_file_name);



%%
big_image=merge_cell_image_1D(roi_list([1:5:15,41]),1,5,0);
imwrite(big_image,[folder,'/roi/combined1'...
            ,'.jpg'],'jpg');
%% and do fine alignment again
if ~exist([folder,'/aligned1'], 'dir')
    mkdir([folder,'/aligned1']);
end

first_frame=Image_coarse_align{frames(1)};
f=figure;
imagesc(first_frame);
axis image;
colormap(gray);

h = imrect;
position = wait(h);
position=floor(position);
cc_rangex=position(2):1:position(2)+position(4);
cc_rangey=position(1):1:position(1)+position(3);
close(f);
shifts2=zeros(length(frames),2);
for i=2:1:length(frames)
    first_frame=Image_coarse_align{frames(i-1)};
    current_frame=Image_coarse_align{frames(i)};
    
    data1=clear_background(first_frame,10,0);
    data2=clear_background(current_frame,10,0);
    data1=data1(cc_rangex,cc_rangey);
    data2=data2(cc_rangex,cc_rangey);
    
    s=dftregistration((fft2(data1)),(fft2(data2)),10);
    shifts2(i,1)=s(3)+shifts2(i-1,1);
    shifts2(i,2)=s(4)+shifts2(i-1,2);
    current_frame=imshift([shifts(i,1) shifts(i,2)],current_frame);
    save_image(current_frame,[folder,'/aligned1/frame'...
            ,num2str(i),'.jpg']);
end

%%
outputVideo = VideoWriter(fullfile(folder,'/shuttle_out'),'MPEG-4');
outputVideo.FrameRate =10;
open(outputVideo)

for i = 1:length(Image_coarse_align)
   writeVideo(outputVideo,abs(Image_coarse_align{i}(1:1088,1:end))/max(Image_coarse_align{i}(:)))
end
close(outputVideo)

%%
outputVideo = VideoWriter(fullfile(folder,'/shuttle_out'),'MPEG-4');
outputVideo.FrameRate = xyloObj.FrameRate;
open(outputVideo)

for i = 1:length(ImageFrames)
   writeVideo(outputVideo,abs(ImageFrames{i})/max(ImageFrames{i}(:)))
end
close(outputVideo)