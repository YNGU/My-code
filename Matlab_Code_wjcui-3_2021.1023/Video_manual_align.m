 % prealignment 从Video_auto_align.m中save the aligned images后开始，手动点位置align图片

%% preshift option
% click on the the characteristic spot
feature_location=zeros(length(frames),2);
pre_shifts=zeros(length(frames),2);
prealignframe=length(frames); %全部图片点出来
prealignframe=79; %需要align多少张图
for k=1:1:prealignframe
    imagesc(ImageFrames{frames(k)});
    axis image
    colormap(inferno);
     title(['frame ' num2str(k)]);
    [x1, y1, but] = ginput(1);
   
    feature_location(k,1)=y1;
    feature_location(k,2)=x1;
    pre_shifts(k,1)=feature_location(1,1)-y1;
    pre_shifts(k,2)=feature_location(1,2)-x1;
end

%% save the prealigned images 最后一步，保存到prealign文件夹，已经prealign好的图片导出
if ~exist([folder,'/prealign'], 'dir')
    mkdir([folder,'/prealign']);
end
shifts=pre_shifts;
x_left =ceil(abs(min(shifts(:,1))))+1;
x_right=ceil(abs(max(shifts(:,1))))+1;
y_top   =ceil(abs(min(shifts(:,2))))+1;
y_bottom=ceil(abs(max(shifts(:,2))))+1;
[x_left x_right y_top y_bottom]
[sizex sizey]=size(ImageFrames{1});
image_till=sizex; % this deals with the time tag
for i=1:1:prealignframe
    large_frame=zeros(x_left+sizex+x_right,y_top+sizey+y_bottom);
    %large_frame(x_left:x_left+image_till-1,y_top:y_top+sizey-1)=...
    %    ImageFrames{frames(i)}(1:image_till,:);
    large_frame(x_left:x_left+image_till-1,y_top:y_top+sizey-1)=...
        ImageFrames{frames(i)}(1:image_till,:);
    Image_coarse_align{i}=imtranslate(large_frame,[shifts(i,2) shifts(i,1)]); %Image_coarse_align是变量
    %     current_frame(x_left:x_left+sizex-1,y_top+image_till:y_top+sizey)=...
    %         ImageFrames{i}(:,image_till:sizey);
    h = fspecial('gaussian', [20 20],1);
    %ImageSum_BL=conv2(Image_coarse_align{i},h,'same');
    ImageSum_BL=Image_coarse_align{i};
    %J = imadjust(Image_coarse_align{i}/255,[0;1], [0;1]);
    J=imwrite_color(ImageSum_BL,1,200);
    imwrite(J,[folder,'/prealign/frame'...
            ,num2str(i),'.jpg'],'jpg');
end

%% align the images to its previous one 给图片align
if ~exist([folder,'/aligned'], 'dir')
    mkdir([folder,'/aligned']);
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

shifts=zeros(length(frames),2);

for i=2:1:20%length(frames)  %改成需要align的图片数量，全部图片是length(frames)
    first_frame=Image_coarse_align{frames(i-1)};
    current_frame=Image_coarse_align{frames(i)};
    
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
%% save the images to compare if the shifts were calculated correctly 对比前后一张图片
if ~exist([folder,'/check_cc'], 'dir')
    mkdir([folder,'/check_cc']);
end

for k=2:1:length(frames)
    Image_shifted=imtranslate(Image_coarse_align{frames(k)},[shifts_adjacent(k,2) shifts_adjacent(k,1)]);
    ImageSum_BL=[Image_coarse_align{frames(k-1)} Image_shifted];
    J = imadjust(ImageSum_BL/255,[0;1], [0;1]);
    imwrite(J,[folder,'/check_cc/frame'...
            ,num2str(k),'.jpg'],'jpg');
end

%% restore shifts 保存偏移量
for i=2:1:length(frames)
    shifts(i,1)=shifts_adjacent(i,1)+shifts(i-1,1);
    shifts(i,2)=shifts_adjacent(i,2)+shifts(i-1,2);
end
%% save the aligned images 最后一步，保存到align文件夹，已经align好的图片导出
x_left =ceil(abs(min(shifts(:,1))))+1;
x_right=ceil(abs(max(shifts(:,1))))+1;
y_top   =ceil(abs(min(shifts(:,2))))+1;
y_bottom=ceil(abs(max(shifts(:,2))))+1;
[x_left x_right y_top y_bottom]
[sizex sizey]=size(Image_coarse_align{1});
image_till=sizex; % this deals with the time tag
for i=1:1:length(frames)
    large_frame=zeros(x_left+sizex+x_right,y_top+sizey+y_bottom);
    %large_frame(x_left:x_left+image_till-1,y_top:y_top+sizey-1)=...
    %    ImageFrames{frames(i)}(1:image_till,:);
    large_frame(x_left:x_left+image_till-1,y_top:y_top+sizey-1)=...
        ImageFrames{frames(i)}(1:image_till,:);
    Image_fine_align{i}=imtranslate(large_frame,[shifts(i,2) shifts(i,1)]);
    %     current_frame(x_left:x_left+sizex-1,y_top+image_till:y_top+sizey)=...
    %         ImageFrames{i}(:,image_till:sizey);
    h = fspecial('gaussian', [20 20],1);
    %ImageSum_BL=conv2(Image_coarse_align{i},h,'same');
    ImageSum_BL=Image_fine_align{i};
    %J = imadjust(Image_coarse_align{i}/255,[0;1], [0;1]);
    J=imwrite_color(ImageSum_BL,1,200);
    imwrite(J,[folder,'/aligned/frame'...
            ,num2str(i),'.jpg'],'jpg');
end
