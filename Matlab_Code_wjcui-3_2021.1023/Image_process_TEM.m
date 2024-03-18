%% 0.clear
clc
clear
%% 1.read image with tif format
imagePath = 'D:\++++++Data\Matlab Data\Exp_Data\Bi liquid\+relat fft and image\791.tif'; % 请替换 'your_image.tif' 为实际图片文件的路径
imageData = imread(imagePath); % 使用imread函数读取tif图片
croppedImageData = imageData;
%% ++.mask particle area using rectangle
% 指定旋转角度（顺时针为正，逆时针为负）
rotationAngle = 13; % 你可以根据需要调整旋转角度

% 对裁剪后的图片进行旋转
rotatedCroppedImageData = imrotate(croppedImageData, rotationAngle, 'bilinear', 'crop');

% 显示旋转后的图片
figure;  
imshow(rotatedCroppedImageData);
title('旋转后的图片');

% 使用imrect工具设置矩形的大小
disp('请手动设置矩形的大小...');
rect = imrect;
wait(rect);

% 获取矩形的位置和大小
position = getPosition(rect);
x1 = round(position(1));
y1 = round(position(2));
width = round(position(3));
height = round(position(4));
x2 = x1 + width - 1;
y2 = y1 + height - 1;

% 创建一个矩形掩模
rectangleMask = zeros(size(rotatedCroppedImageData));
rectangleMask(y1:y2, x1:x2) = 1;

% 对旋转后的裁剪图片应用矩形掩模
maskedRegion = double(rotatedCroppedImageData) .* rectangleMask;

% 显示应用了矩形掩模的图像
figure;
imshow(uint8(maskedRegion));
title('应用了矩形掩模的图像');
%% ++.hanning window fft for mask area
% 将应用了圆形掩模的图像应用Hanning窗
hanningWindow = hann(size(maskedRegion, 1)) * hann(size(maskedRegion, 2))';
maskedRegion = maskedRegion .* hanningWindow;

% 进行FFT变换
fftResult = fft2(maskedRegion);

% 显示FFT结果的幅度谱
figure;
imshow(log(abs(fftshift(fftResult)) + 1), []);
title('圆形掩模FFT变换的幅度谱');
%% 

