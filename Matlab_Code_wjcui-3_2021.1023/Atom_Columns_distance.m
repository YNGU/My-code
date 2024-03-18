%% clear all
clear
clc
%% read the dm3 image
dm3_file_name='D:\++++++Data\Matlab Data\Exp_Data\LaFeSi_sunyue\1607 8.1 Mx HAADF Nano STEM 20220605 Diffraction.dm3';
ImageSum_R=ReadDMFile(dm3_file_name); 
ImageSum_R=double(ImageSum_R);
%% 1.read tif image 
ImageSum_R=imread('D:\++++++Data\Matlab Data\Exp_Data\BiFeO\1519 8.1 Mx HAADF Nano Diffraction STEM 20220312.tif');
ImageSum_R(:,:,4)=[]; %四通道必要时再用
ImageSum_R=rgb2gray(ImageSum_R);
ImageSum_R=double(ImageSum_R);
%% 2.look image
figure;
ImageSum_R_blurred=imgaussfilt(ImageSum_R,1);

imagesc(ImageSum_R_blurred);
axis image;
colormap(gray);
%% 3.crop the image if necessary
figure;
imagesc(ImageSum_R);
axis image;
colormap(hot);

h = imrect;

position = wait(h);
position=floor(position);
image_rangex=position(2):1:position(2)+position(4);
image_rangey=position(1):1:position(1)+position(3);
ImageSum_R=ImageSum_R(image_rangex,image_rangey);
imagesc(ImageSum_R);
axis image
%% 4.gaussian distribution
[fitresult_N,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,object_index,mass_center_N]=...
    find_atomic_columns(ImageSum_R,0,0.1,6,1,1,[0 350],0,0,1);%ImageSum_R,0(变值为0-6),0.1（值越大，阈值越大，原子分的越精细）,6,1,1,[100 400] %区间设定两个峰为取亮暗原子，前峰为暗原子，后峰亮原子，取想要的原子区间；%find-atomic-columns: template=sign*fspecial('gaussian',template_size*2-1,sigma)里面template-size数值可以更改小点，找到更弱的点但是噪音也可能被标记出来。
%% 5.find the atom column locations using normalized cross correlation data
% function [fitresult,oimage,zfit, fiterr, zerr, resnorm,...
% rr,image1,image2,object_index,mass_center,C,D,StartPoint,h_area]=...
% find_atomic_columns(raw_image,sigma,threshold,max_peak_num,sign,style,...
% area_threshold,initial_values,fit_shift,verbose);

% the main output 'fitresult_N' has a cell structure,
% each cell is a 1x7 array [amp, ang, sx, sy, xo, yo, zo] containing the peak fitting result
% amp: amplitude
% ang: rotation angle of the two main axes
% sx: sigma along the first main axis
% sy: sigma along the second main axis
% xo: x coordinate of the peak center
% yo: y coordinate of the peak center
% zo: background intensity

% Meaning of inputs of 'find_atom_clolumns'
% raw_image: The input STEM image
% sigma: sigma of the gaussian distribution for normalized cross-correlation,
% when the number is 0, the program decides the best sigma for ncc
% threshold: threshold to separate the atom columns
% max_peak_num: the limit of peak numbers
% sign: 1: find the peaks; 2: find the valleys (for ABF)
% style: 1: use ncc data to fit the peak; 2: use experimental data
% area_threshold: only areas within the area_threshold range are used for fitting 
% initial_values: starting values for peak fitting, can be set as 0
% fit_shift: for future use, 0
% verbose: show the ncc threshold map and area size histogram
[fitresult_N,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,object_index,mass_center_N]=...
    find_atomic_columns(ImageSum_R,0,0.1,6000,1,1,[0 250],0,0,1);%ImageSum_R,0,0.1与上个语句相同，6000为计算点的数目
%% 5.run this instead
fitresult_E=fitresult_N;

%% 6.distance histogram calculated from the fitting result
figure;
[xydist_E,h_E]=position_analysis(fitresult_E,ImageSum_R,300,30);
plot(h_E);

%% 7.calculate PSD for the example image
figure;
[Iproj,Iavg,Istd]=project_image_RD(ImageSum_R,100,-90:1:90);
plot(Istd);

%% 8.we can then use the point-PSD to find exactly locations of peaks at roughly 85 degrees and -5 degrees (175 degrees in the plot) 
figure;
[d_E,proj_acc_E,projx_E,peak_index_E,row_map_E,col_map_E,mini_E,row_stat_E,col_stat_E,coord_angle_E]...
    =assign_xy_to_peaks(ImageSum_R,fitresult_E,200,1,91,1,1); 
plot(projx_E);

%% 9.note the output from matlab says 'find image aligned a long xxx degree at index yyy', use the index yyy to plot the projected profile

figure;
subplot(2,1,1); plot(proj_acc_E(45,:));
subplot(2,1,2); plot(proj_acc_E(145,:));

%% 10.the matrix representation is stored in 'mini_E'
% mini_E is a 2D matrix with each node containing the index of the peak
% fitting result in fitresult_R

figure;
imagesc(mini_E>0);
%daspect([1 sqrt(2) 1]);
daspect([1 1 1]);
colormap(gray);

%% 11.calculate the intensity map easily from the matrix representation
col_int_E=get_col_int( ImageSum_R,fitresult_E,0,0);%col_int是原子颜色强度的变量
% the circle around each atom coloumn shows the intensity
figure;
imagesc(ImageSum_R);
axis image
axis off
colormap(gray);
hold on;

[sx,sy]=size(mini_E);
hold all;
jets=cool(256);%jet是colormap的一种，从蓝色到红色的标准标尺,可以变jet为inferno，hot，cool，hsv等colormap
upper=max(col_int_E);%标尺中颜色的最高值
lower=min(col_int_E);%标尺中颜色的最低值
upper=1;%可以改变最高最低值来尝试图像的效果
lower=0.8;
%lower=1 %只标亮原子的时候用
for i=1:1:sx
    for j=1:1:sy
        if mini_E(i,j) == 0
            continue;
        end
        p1=mini_E(i,j);
        if col_int_E(p1)==0
            continue;
        end
        color_temp=round((col_int_E(p1)-lower)/(upper-lower)*256);%color的算法，可以变此公式把颜色区间变大
        if(color_temp<1) color_temp=1; end
        if(color_temp>256) color_temp=256; end
        plot(fitresult_E{p1}(5), fitresult_E{p1}(6),'or','color',jets(color_temp,:),'markersize',6,'linewidth',8);%变linewidth为实心圆还是圈
    end
end
title('intensity map','fontsize',20);
%% 11-save as
imwrite(col_int_E,'E:\++++++Matlab-Data\Exp data\Hu Zhi-Yi\STEM HAADF 20200611 1404-GF-intensity.tif');
%% 12.separate A site and B site
[mini_s_E,mini_w_E,mini_s_int_E,mini_w_int_E,index_s_E,index_w_E] =...
    mini_separate_weak_strong(mini_E,fitresult_E,ImageSum_R,peak_index_E,col_int_E,[2 2]);
%% 13.
figure;
subplot(1,2,1), imagesc(mini_s_E);axis image
subplot(1,2,2), imagesc(mini_w_E);axis image

%% 14.distance of weak atoms
f=figure;
scale_bar=1;
f.PaperPosition=[1 1 7 7];
lattice_x=draw_lines_from_mini_peak_pair(mini_w_E(1:end,1:end),fitresult_E,...
    ImageSum_R,scale_bar,[0 2],2,4,[0 0],[23.5 25.3]);
axis off
%print(f,'-dpng','-r600',['lattice_x.png']);

%% distance of strong atoms
f=figure;
scale_bar=1;
f.PaperPosition=[1 1 7 7];
lattice_x=draw_lines_from_mini_peak_pair(mini_s_E(1:end,1:end),fitresult_E,...
    ImageSum_R,scale_bar,[0 2],2,4,[0 0],[74 74.1]);%区间按照ans输出的值变化
axis off
%print(f,'-dpng','-r600',['lattice_x.png']);
%% 15.average distance of p14
lattice1=lattice_x';
s1=size(lattice1);
avg_c=zeros(s1(1),1);
std_c=zeros(s1(1),1);
for i=1:1:s1(1)
    temp=lattice1(i,:);
    temp=temp(temp>0);
    avg_c(i)=mean(temp);
    std_c(i)=std(temp);
end
%plot(avg_c,'o');
errorbar(avg_c,std_c,'or');
%% save mat file
mat_file_name='E:\Matlab Data\Experiment-Data\STEM HAADF 20190915 0028.mat'%把文件名后缀改成mat,要打开该文件时load（‘文件名’）就可以
save(mat_file_name)