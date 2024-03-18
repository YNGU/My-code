clc
clear
%% 1. read image  读文件
%ImageSum_R=imread('E:\GeTe\2020 10 21 GeSbTe HDAAF\IFFT of Untitled 2.tif');
%ImageSum_R=imread('E:\GeTe\20201022-GeSbTe-tif\STEM HAADF 1706 20201021 8.1 Mx B DCFI(HAADF).tif');
%ImageSum_R=imread('E:\GeTe\20201022-GeSbTe-tif\STEM HAADF 1722 20201021 8.1 Mx B DCFI(HAADF).tif');
%ImageSum_R=imread('E:\GeTe\20201022-GeSbTe-tif\STEM HAADF 1737 20201021 5.7 Mx B DCFI(HAADF).tif');
%ImageSum_R=imread('E:\GeTe\20201022-GeSbTe-tif\1737.tif');
%ImageSum_R=imread('E:\GeTe\20201022-GeSbTe-tif\STEM HAADF 1736 20201021 5.7 Mx B DCFI(HAADF).png');
%ImageSum_R=imread('E:\GeTe\20201022-GeSbTe-tif\1735.tif');
%ImageSum_R=imread('E:\GeTe\20201022-GeSbTe-tif\STEM HAADF 1751 20201021 5.7 Mx B DCFI(HAADF).tif');
%ImageSum_R=imread('E:\sample\1.tif');
%ImageSum_R=imread('E:\sample\yuyimeng\4cpg.tif');
%ImageSum_R=imread('E:\sample\mengxiangyu\STEM HAADF 1719 20201021 4.0 Mx HAADFc.tif');
%ImageSum_R=imread('E:\sample\mengxiangyu\GeTe\STEM HAADF 0620 20201221 4.0 Mx c.tif');
%ImageSum_R=imread('E:\sample\mengxiangyu\GeTe\STEM HAADF 0636 20201221 5.7 Mx.tif');
%ImageSum_R=imread('E:\sample\mengxiangyu\GeTe\STEM HAADF 0637 20201221 5.7 Mx.tif');
%ImageSum_R=imread('E:\sample\mengxiangyu\GeTe1\4.tif');
%ImageSum_R=imread('G:\polarization\20220730 A2\1\0923 8.1 Mx HAADF Nano STEM 20220730 Diffraction c.tif');
ImageSum_R=imread('D:\++++++Data\Matlab Data\Exp_Data\test\2251 5.7 Mx HAADF STEM.tif');
%ImageSum_R=imread('D:\++++++Data\Matlab Data\Exp_Data\Test\0856 5.7 Mx HAADF Nano STEM 20220730 Diffraction 0001 c.tif');
%% 2. change unit to double 转格式
ImageSum_R(:,:,4)=[];
%ImageSum_R=uint8(ImageSum_R);
ImageSum_R=rgb2gray(ImageSum_R);
ImageSum_R=double(ImageSum_R);
%% 3. look image 看图像
figure;
ImageSum_R_blurred=imgaussfilt(ImageSum_R,1);
imagesc(ImageSum_R_blurred);
axis image;
colormap(gray);
%% crop the image if necessary 可以不裁剪
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
%% 4. gaussian distribution  强度分布
[fitresult_N,oimage,zfit, fiterr, zerr, resnorm, rr,image1,image2,object_index,mass_center_N]=...
find_atomic_columns(ImageSum_R,0,0.1,6,1,1,[0,200],0,0,1);
%(raw_image,sigma,threshold,max_peak_num,sign,style,area_threshold,initial_values,fit_shift,verbose)
%% 5. find the atom column locations using normalized cross correlation data 找点数
% function [fitresult,oimage,zfit, fiterr, zerr, resnorm,...
% rr,image1,image2,object_index,mass_center,C,D,StartPoint,h_area]=...
% find_atomic_columns(raw_image,sigma,threshold,max_peak_num,sign,style,...
% area_threshold,initial_values,fit_shift,verbose);

% the main output 'fitresult_N' has a cell structure,
% each cell is a 1x7 array [amp, ang, sx, sy, xo, yo, zo] containing the peak fitting result
% amp: amplitude%振幅
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
    find_atomic_columns(ImageSum_R,0,0.1,20000,1,1,[0 200],0,0,1);

%% 6. find the atom column locations using experimental RevSTEM data skip it 找点，这一步后面没用到
% style is set to 2 to use experimental data for fitting
% starting points are set to be fitresult_N
[fitresult_E,oimage_E,zfit_E, fiterr_E, zerr_E, resnorm_E, rr_E,image1_E,image2_E,startpoint_E,mass_center_E]=...
    find_atomic_columns(ImageSum_R,0,0.1,20000,1,1,[0 350],fitresult_N,0,0);


%% 存成mat文件
%mat_file_name='E:\GeTe\20201022-GeSbTe-tif\STEM HAADF 1706 20201021 8.1 Mx B DCFI(HAADF).mat';
%mat_file_name='E:\GeTe\20201022-GeSbTe-tif\STEM HAADF 1706 20201021 8.1 Mx B DCFI(HAADF).mat';
%mat_file_name='E:\GeTe\20201022-GeSbTe-tif\1737.mat';
%mat_file_name='E:\GeTe\20201022-GeSbTe-tif\STEM HAADF 1736 20201021 5.7 Mx B DCFI(HAADF).mat';
%mat_file_name='E:\GeTe\20201022-GeSbTe-tif\1737 2.mat';
%mat_file_name='E:\GeTe\20201022-GeSbTe-tif\1735.mat';
%mat_file_name='E:\sample\yuyimeng\4cpg.mat';
%mat_file_name='E:\sample\mengxiangyu\GeTe\STEM HAADF 0636 20201221 5.7 Mx.tif';
%mat_file_name='D:\++++++Data\Matlab Data\Exp_Data\Test\0856 5.7 Mx HAADF Nano STEM 20220730 Diffraction 0001 c.mat';
%% 存文件
save(mat_file_name);
%% 7. distance histogram calculated from the fitting result 
figure;
[xydist_E,h_E]=position_analysis(fitresult_E,ImageSum_R,300,30);
plot(h_E);

%% 8. Peak_position  看点是否找的有问题
d_f=length(fitresult_E);
Peak_position=zeros(2,d_f);

for i=1:d_f
    Peak_position(1,i)=fitresult_E{i}(5);
    Peak_position(2,i)=fitresult_E{i}(6);
   
end
figure;
imagesc(ImageSum_R);
axis image
axis off
colormap(gray);
hold on;
for i=1:d_f
 plot(Peak_position(1,i),Peak_position(2,i),'.','color','red','markersize',10,'linewidth',2);
end
%% 9. calculate PSD for the example image 找投影矢量度数
figure;
[Iproj,Iavg,Istd]=project_image_RD(ImageSum_R,100,-90:1:90); %不用改
plot(Istd);

%% 10. we can then use the point-PSD to find exactly locations of peaks at roughly 85 degrees and -5 degrees (175 degrees in the plot) 
figure;
[d_E,proj_acc_E,projx_E,peak_index_E,row_map_E,col_map_E,mini_E,row_stat_E,col_stat_E,coord_angle_E]...
    =assign_xy_to_peaks(ImageSum_R,fitresult_E,200,2,90,3,1);  % 换1，91为上面两个90度峰的值
plot(projx_E);                                  

%% 11. note the output from matlab says 'find image aligned a long xxx degree at index yyy', use the index yyy to plot the projected profile

figure;
subplot(2,1,1); plot(proj_acc_E(23,:));
subplot(2,1,2); plot(proj_acc_E(137,:));

%% 12. the matrix representation is stored in 'mini_E'
% mini_E is a 2D matrix with each node containing the index of the peak
% fitting result in fitresult_R

figure;
imagesc(mini_E>0);
%daspect([1 sqrt(2) 1]);
daspect([1 1 1]);
colormap(gray);
%% 13.
col_int_E=get_col_int( ImageSum_R,fitresult_E,0,0);

%% 14. separate use intensityseparate use intensity
psize=length(fitresult_E);
[ImageX, ImageY]=size(ImageSum_R);
[sx,sy]=size(mini_E);
mini_int=zeros(sx,sy);
mini_s=zeros(sx,sy);
mini_w=zeros(sx,sy);
mini_s_int=zeros(sx,sy);
mini_w_int=zeros(sx,sy);

%for i=1:1:psize
 %   if(peak_index(i,1)>0 && peak_index(i,2)>0)
 %       mini_int(peak_index(i,1),peak_index(i,2))=col_int(i);
 %   end
%end
for i=1:1:sx
    for j=1:1:sy
        if mini_E(i,j)==0
            continue;
        end
         p2=mini_E(i,j);
         a1=fitresult_E{p2}(1)+fitresult_E{p2}(7);
         mini_int(i,j)=a1;
        end
    end
%% 14. 
hist(col_int_E); %统计亮原子与暗原子的强度位置，找到两者之间的阈值
%% 15. 
[sx,sy]=size(mini_E);
mini_s_E=zeros(sx,sy);
mini_s_int=zeros(sx,sy);
for i=1:1:sx
    for j=1:1:sy
        if mini_int(i,j)>0.75 % 找亮原子A位置，找到上一步的阈值，输入该位置
            mini_s_E(i,j)=mini_E(i,j);
            mini_s_int(i,j)=mini_int(i,j);
        end
    end
end

mini_w_E=zeros(sx,sy);
mini_w_int=zeros(sx,sy);
for i=1:1:sx
    for j=1:1:sy
        if mini_int(i,j)<0.75 % 找暗原子B的位置
            mini_w_E(i,j)=mini_E(i,j);
            mini_w_int(i,j)=mini_int(i,j);
        end
    end
end
%% 15. quiver map 得到偏移map
figure;
imagesc(ImageSum_R);
axis image
axis off
colormap(gray);
hold on;


[sx,sy]=size(mini_E);
hold all;
dx=zeros(sx,sy);
dy=zeros(sx,sy);
d_abs=zeros(sx,sy);
quiver_x=[];
quiver_y=[];
quiver_dx=[];
quiver_dy=[];
count=1;
for i=2:1:sx-1
    for j=2:1:sy-1
        if mini_w_E(i,j) == 0
            continue;
        end
        p1=mini_w_E(i,j);
        if mini_s_E(i-1,j+1)==0 || mini_s_E(i+1,j-1)==0 || mini_s_E(i+1,j+1)==0 || mini_s_E(i-1,j-1)==0 %根据图像的（i，j）排列规则改变参数，图像不正改参数
            continue;
        end
        
        ps_1=mini_s_E(i-1,j+1); %图像如果不正需要改参数，i和j
        ps_2=mini_s_E(i+1,j-1);
        ps_3= mini_s_E(i+1,j+1);
        ps_4=mini_s_E(i-1,j-1);
        deltax=fitresult_E{p1}(5)-(fitresult_E{ps_1}(5)+fitresult_E{ps_2}(5)+fitresult_E{ps_3}(5)+fitresult_E{ps_4}(5))/4;
        deltay=fitresult_E{p1}(6)-(fitresult_E{ps_1}(6)+fitresult_E{ps_2}(6)+fitresult_E{ps_3}(6)+fitresult_E{ps_4}(6))/4;
        dx(i,j)=deltax;
        dy(i,j)=deltay;
        d_abs(i,j)=sqrt(deltax^2+deltay^2);
        quiver_x(count)=fitresult_E{p1}(5);
        quiver_y(count)=fitresult_E{p1}(6);
        quiver_dx(count)=deltax;
        quiver_dy(count)=deltay;
        count=count+1;
    end
end
quiver(quiver_x,quiver_y,15*quiver_dx,15*quiver_dy,1.4,'LineWidth',1.4,'color','y'); %25*quiver_dx和dy改变前面乘值改箭头长度大小
title('quiver','fontsize',20);
%% a strain map
[strain_a,maxtri_a]=a_strain(fitresult_E,mini_w_E,mini_s_E,ImageSum_R,1,1);
figure;
imagesc(ImageSum_R);                                                                                        
axis image
axis off
colormap(gray);
hold on;


[sx,sy]=size(mini_E);
hold all;
jets=mycolor;
%de_s=strain_a(strain_a>0);
%upper=max(strain_a);
%lower=min(strain_a);
strain_a_E=strain_a(strain_a>0);
upper=max(strain_a_E);
lower=min(strain_a_E);
for i=2:1:sx-1
    for j=2:1:sy-1
        if mini_w_E(i,j) == 0
            continue;
        end
        if mini_s_E(i-1,j)==0 || mini_s_E(i+1,j)==0 || mini_s_E(i,j+1)==0 || mini_s_E(i,j-1)==0
            continue;
        end
        if maxtri_a(i,j)==0
            continue;
        end
        p1=mini_w_E(i,j);
        %if strain_a(p1)==0
           % continue;
        %end
      
        color_temp=round((strain_a(p1)-lower)/(upper-lower)*64);
        if(color_temp<1) color_temp=1; end
        if(color_temp>256) color_temp=64; end
        plot(fitresult_E{p1}(5), fitresult_E{p1}(6),'.','color',jets(color_temp,:),'markersize',35,'linewidth',2);
    end
end
title('a strain map','fontsize',20);
hold off
%% colorbar
figure;
axis image
axis off
colormap(mycolor);
colorbar('Ticks',[0,1],'TickLabels',{'1.40','1.47'});
%% c strain map
[strain_c,maxtri_c]=c_strain(fitresult_E,mini_w_E,mini_s_E,ImageSum_R,1,1);
figure;
imagesc(ImageSum_R);
axis image
axis off
colormap(gray);
hold on;


[sx,sy]=size(mini_E);
hold all;
jets=mycolor;
strain_c_E=strain_c(strain_c>0);
upper=max(strain_c_E);
lower=min(strain_c_E);

for i=2:1:sx-1
    for j=2:1:sy-1
        if mini_w_E(i,j) == 0
            continue;
        end
        if mini_s_E(i-1,j)==0 || mini_s_E(i+1,j)==0 || mini_s_E(i,j+1)==0 || mini_s_E(i,j-1)==0
            continue;
        end
        if maxtri_c(i,j)==0
            continue;
        end
        p1=mini_w_E(i,j);
        %if strain_c_E(p1)==0
            %continue;
        %end
        
        color_temp=round((strain_c(p1)-lower)/(upper-lower)*64);
        if(color_temp<1) color_temp=1; end
        if(color_temp>64) color_temp=64; end
        plot(fitresult_E{p1}(5), fitresult_E{p1}(6),'.','color',jets(color_temp,:),'markersize',35,'linewidth',2);
    end
end
title('c strain map','fontsize',20);
hold off
%% colorbar
figure;
axis image
axis off
colormap(mycolor);
colorbar('Ticks',[0,1],'TickLabels',{'1.30','1.38'});
%%
%maxtri_a_R=
cir1=generate2Dcircle(5,5,10);
ima_angle_blur=conv2(maxtri_a,cir1,'same');
showarealimage(maxtri_a);
colormap(mycolor)
colorbar

%% calculate the intensity map of strong atoms easily from the matrix represe
col_int_E=get_col_int( ImageSum_R,fitresult_E,0,0);
% the circle around each atom coloumn shows the intensity
figure;
imagesc(ImageSum_R);
axis image
axis off
colormap(gray);
hold on;

[sx,sy]=size(mini_E);
hold all;
%jets=winter(200);
jets=mycolor;
for i=1:1:sx
    for j=1:1:sy
        if mini_s_E(i,j) == 0
            continue;
        end
        p1=mini_s_E(i,j);
      
        int_s=mini_s_int(mini_s_int>0);
        upper=max(int_s);
        lower=min(int_s); 
        color_temp=round((col_int_E(p1)-lower)/(upper-lower)*64);
        if(color_temp<1), color_temp=1; end
        if(color_temp>64), color_temp=64; end
        plot(fitresult_E{p1}(5), fitresult_E{p1}(6),'.','color',jets(color_temp,:),'markersize',35,'linewidth',2);
    end
end
%% colorbar
figure;
axis image
axis off
colormap(mycolor);
colorbar('Ticks',[0,1],'TickLabels',{'180','270'});
%% calculate the intensity map of weak atoms easily from the matrix represe
%col_int_E=get_col_int( ImageSum_R,fitresult_E,0,0);
% the circle around each atom coloumn shows the intensity
figure;
imagesc(ImageSum_R);
axis image
axis off
colormap(gray);
hold on;

[sx,sy]=size(mini_E);
hold all;
jets=mycolor;

for i=1:1:sx
    for j=1:1:sy
        if mini_w_E(i,j) == 0
            continue;
        end
        p1=mini_w_E(i,j);
        
        int_w=mini_w_int(mini_w_int>0);
        %upper=max(int_w);
        upper=150;
        lower=min(int_w);
        color_temp=round((col_int_E(p1)-lower)/(upper-lower)*64);
        if(color_temp<1), color_temp=1; end
        if(color_temp>64), color_temp=64; end
        plot(fitresult_E{p1}(5), fitresult_E{p1}(6),'.','color',jets(color_temp,:),'markersize',35,'linewidth',2);
    end
end
title('intensity map','fontsize',20);
hold off

%% colorbar
figure;
axis image
axis off
colormap(mycolor);
colorbar('Ticks',[0,1],'TickLabels',{'95','150'});
%%
%[exx,exy,eyx,eyy,u1x,u1y,u2x,u2y]=strain_analysis_mini(ImageSum_R,fitresult_E,mini_s_E,[2 0],[0 2],4);
figure;
imagesc(exx);
axis image
axis off
colormap(mycolor);
figure;
imagesc(exy);
axis image
axis off
colormap(mycolor);
figure;
imagesc(eyx);
axis image
axis off
colormap(mycolor);
figure;
imagesc(eyy);
axis image
axis off
colormap(mycolor);
%%
psize=length(fitresult_E);
[ImageX, ImageY]=size(ImageSum_R);
[sx,sy]=size(mini_E);
a1=[0 2];
a2=[2 0];
exx=zeros(sx,sy);
exy=zeros(sx,sy);
eyx=zeros(sx,sy);
eyy=zeros(sx,sy);

u1x=zeros(sx,sy);
u1y=zeros(sx,sy);
u2x=zeros(sx,sy);
u2y=zeros(sx,sy);
for i=1:1:sx
    for j=1:1:sy
        if(mini_E(i,j)==0)
            continue;
        end
        p1=mini_E(i,j);
        x1=fitresult_E{p1}(6);
        y1=fitresult_E{p1}(5);
        
        k=i+a1(1);
        l=j+a1(2);
        if (k<1 || l<1 || k>sx || l>sy)
            continue;
        end
        if(mini_E(k,l)==0)
            continue;
        end
        p2=mini_E(k,l);
        x2=fitresult_E{p2}(6);
        y2=fitresult_E{p2}(5);
        u1x(i,j)=x2-x1;
        u1y(i,j)=y2-y1;
        
    end
end
%%
for i=1:1:sx
    for j=1:1:sy
        if(mini_E(i,j)==0)
            continue;
        end
        p1=mini_E(i,j);
        x1=fitresult_E{p1}(6);
        y1=fitresult_E{p1}(5);
        
        k=i+a2(1);
        l=j+a2(2);
        if (k<1 || l<1 || k>sx || l>sy)
            continue;
        end
        if(mini_E(k,l)==0)
            continue;
        end
        p2=mini_E(k,l);
        x2=fitresult_E{p2}(6);
        y2=fitresult_E{p2}(5);
        u2x(i,j)=x2-x1;
        u2y(i,j)=y2-y1;
        
    end
end
%%
u1x_avg=0;
u1y_avg=0;
count=0;
for i=1:1:sx
    for j=1:1:sy
        if(mini_E(i,j)==0)
           continue;
        end
        if u1x(i,j)==0 && u1y(i,j)==0
            continue;
        end
        u1x_avg=u1x_avg+u1x(i,j);
        u1y_avg=u1y_avg+u1y(i,j);
        count=count+1;
    end
end
u1x_avg=u1x_avg/count;
u1y_avg=u1y_avg/count;

u2x_avg=0;
u2y_avg=0;
count=0;
for i=1:1:sx-2
    for j=1:1:sy-2
        if(mini_E(i,j)==0)
            continue;
        end
        if u2x(i,j)==0 && u2y(i,j)==0
            continue;
        end
        u2x_avg=u2x_avg+u2x(i,j);
        u2y_avg=u2y_avg+u2y(i,j);
        count=count+1;
    end
end
u2x_avg=u2x_avg/count;
u2y_avg=u2y_avg/count;

[u1x_avg u1y_avg u2x_avg u2y_avg]

A=[u1x_avg u2x_avg; u1y_avg u2y_avg];
G=inv(A);

a1x=u1x*G(1,1)+u1y*G(1,2);
a1y=u1x*G(2,1)+u1y*G(2,2);

a2x=u2x*G(1,1)+u2y*G(1,2);
a2y=u2x*G(2,1)+u2y*G(2,2);

exx=a1x;
exy=a1y;
eyx=a2x;
eyy=a2y;
%%
for i=1:1:sx
    for j=1:1:sy
        if(exx(i,j)>0)
            exx(i,j)=exx(i,j)-1;
        end
        if(eyy(i,j)>0)
            eyy(i,j)=eyy(i,j)-1;
        end
        
    end
end
%%
%angles(i,j)=acos(abs(deltay)/d_abs(i,j))*180/pi;
figure;
imagesc(ImageSum_R);
axis image
axis off
colormap(gray);
hold on;

[sx,sy]=size(mini_E);
hold all;
jets=mycolor;
dx=zeros(sx,sy);
dy=zeros(sx,sy);
d_abs=zeros(sx,sy);
quiver_x=[];
quiver_y=[];
quiver_dx=[];
quiver_dy=[];
count=1;
for i=2:1:sx-1
    for j=2:1:sy-1
        if mini_w_E(i,j) == 0
            continue;
        end
        p1=mini_w_E(i,j);
        if mini_s_E(i-1,j)==0 || mini_s_E(i+1,j)==0 || mini_s_E(i,j+1)==0 || mini_s_E(i,j-1)==0
            continue;
        end
        
        ps_1=mini_s_E(i-1,j);
        ps_2=mini_s_E(i+1,j);
        ps_3= mini_s_E(i,j+1);
        ps_4=mini_s_E(i,j-1);
        deltax=fitresult_E{p1}(5)-(fitresult_E{ps_1}(5)+fitresult_E{ps_2}(5)+fitresult_E{ps_3}(5)+fitresult_E{ps_4}(5))/4;
        deltay=fitresult_E{p1}(6)-(fitresult_E{ps_1}(6)+fitresult_E{ps_2}(6)+fitresult_E{ps_3}(6)+fitresult_E{ps_4}(6))/4;
        dx(i,j)=deltax;
        dy(i,j)=deltay;
        d_abs(i,j)=sqrt(deltax^2+deltay^2);
        angles(i,j)=acos(deltay/d_abs(i,j))*180/pi;
        quiver_x(count)=fitresult_E{p1}(5);
        quiver_y(count)=fitresult_E{p1}(6);
        quiver_dx(count)=deltax;
        quiver_dy(count)=deltay;
        count=count+1;
    end
end
for i=2:1:sx-1
    for j=2:1:sy-1
        if mini_w_E(i,j) == 0
            continue;
        end
        p1=mini_w_E(i,j);
        if mini_s_E(i-1,j)==0 || mini_s_E(i+1,j)==0 || mini_s_E(i,j+1)==0 || mini_s_E(i,j-1)==0
            continue;
        end
        %if angles(i,j)==0
           % continue;
       % end
        %upper=180;
        %lower=-180;
        d_abs_nz=d_abs(d_abs>0); % remove data points that are zero
        upper=max(d_abs_nz);
        lower=min(d_abs_nz);
        color_temp=round((d_abs(i,j)-lower)/(upper-lower)*64);
        if(color_temp<1) color_temp=1; end
        if(color_temp>64) color_temp=64; end
        plot(fitresult_E{p1}(5), fitresult_E{p1}(6),'.','color',jets(color_temp,:),'markersize',35,'linewidth',2);
    end
end
%% colorbar
figure;
axis image
axis off
colormap(mycolor);
colorbar('Ticks',[0,1],'TickLabels',{'0.06','1.13'});
%% quiver intensity map
figure;
imagesc(ImageSum_R);
axis image
axis off
colormap(gray);
hold on;

[sx,sy]=size(mini_E);
hold all;
dx=zeros(sx,sy);
dy=zeros(sx,sy);
d_abs=zeros(sx,sy);
quiver_x=[];
quiver_y=[];
quiver_dx=[];
quiver_dy=[];
count=1;
for i=2:1:sx-1
    for j=2:1:sy-1
        if mini_w_E(i,j) == 0
            continue;
        end
        p1=mini_w_E(i,j);
        if mini_s_E(i-1,j)==0 || mini_s_E(i+1,j)==0 || mini_s_E(i,j+1)==0 || mini_s_E(i,j-1)==0
            continue;
        end
        
        ps_1=mini_s_E(i-1,j);
        ps_2=mini_s_E(i+1,j);
        ps_3= mini_s_E(i,j+1);
        ps_4=mini_s_E(i,j-1);
        deltax=fitresult_E{p1}(5)-(fitresult_E{ps_1}(5)+fitresult_E{ps_2}(5)+fitresult_E{ps_3}(5)+fitresult_E{ps_4}(5))/4;
        deltay=fitresult_E{p1}(6)-(fitresult_E{ps_1}(6)+fitresult_E{ps_2}(6)+fitresult_E{ps_3}(6)+fitresult_E{ps_4}(6))/4;
        dx(i,j)=deltax;
        dy(i,j)=deltay;
        d_abs(i,j)=sqrt(deltax^2+deltay^2);
        quiver_x(count)=fitresult_E{p1}(5);
        quiver_y(count)=fitresult_E{p1}(6);
        quiver_dx(count)=deltax;
        quiver_dy(count)=deltay;
        count=count+1;
    end
end
for i=2:1:sx-1
    for j=2:1:sy-1
        if mini_w_E(i,j) == 0
            continue;
        end
        p1=mini_w_E(i,j);
        if mini_s_E(i-1,j)==0 || mini_s_E(i+1,j)==0 || mini_s_E(i,j+1)==0 || mini_s_E(i,j-1)==0
            continue;
        end
        
        d_abs_nz=d_abs(d_abs>0); % remove data points that are zero
        upper=max(d_abs_nz);
        lower=min(d_abs_nz);
        color_temp=round(d_abs_nz(i,j)-lower/(upper-lower)*64);
        if(color_temp<1) color_temp=1; end
        if(color_temp>64) color_temp=64; end
        plot(fitresult_E{p1}(5), fitresult_E{p1}(6),'.','color',jets(color_temp,:),'markersize',45,'linewidth',2);
    end
end
