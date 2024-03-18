%% 0.clear
clc
clear
load('G:\++++++data_process\+bigNW_2241_piexl size 0.034nm\output_frames\output_dis_angle.mat')
%% 1.read image with tif format 输入tif图片
imagePath = 'G:\++++++data_process\+bigNW_2241_piexl size 0.034nm\output_frames\Aligned 2241 460.0 kx Ceta Camera0553.tif'; % 请替换 'your_image.tif' 为实际图片文件的路径
imageData = imread(imagePath); % 使用imread函数读取tif图片
%% 2.crop the import image 大致裁剪一下图片
imageData = imread(imagePath); % 使用imread函数读取tif图片
imshow(imageData);
title('导入的图片'); % 显示导入的图片
croppedImage = imcrop; % 使用imcrop函数手动裁剪图像
figure;  
imshow(croppedImage);
title('裁剪后的图片'); % 显示裁剪后的图像
close
%% 3.rotate the image 旋转图片，表面转水平，再次裁剪。要设旋转角度！
% 指定图片旋转的角度参数
rotationAngle = 74;
% 使用imrotate函数旋转图片
rotatedImage = imrotate(croppedImage, rotationAngle, 'bilinear', 'crop');
% 显示旋转后的图片
imshow(rotatedImage);
title('旋转后的图片');
% 启动交互式裁剪工具
h = imrect;
% 等待用户手动选择裁剪区域并获取裁剪框的位置信息
wait(h);
cropPosition = getPosition(h);
% 使用imcrop函数根据用户选择的位置裁剪旋转后的图片
finalManuallyCroppedImage = imcrop(rotatedImage, cropPosition);
% 显示最终手动裁剪后的图片
imshow(finalManuallyCroppedImage);
title('最终手动裁剪后的图片');
close
%% 4.save the final corpped image 保存最终裁剪的图片
% 获取上一级目录路径
parentFolder = fileparts(imagePath);
% 获取原始文件名的后四位
[~, filename, ~] = fileparts(imagePath);
suffix = filename(end-3:end);
% 创建新文件夹，文件夹名为原始文件名的后四位
newFolderName = suffix;
newFolderPath = fullfile(parentFolder, newFolderName);
% 如果文件夹不存在，则创建
if ~exist(newFolderPath, 'dir')
    mkdir(newFolderPath);
end
% 保存裁剪后的图像到新文件夹
imwrite(finalManuallyCroppedImage, fullfile(newFolderPath, 'cropped_image.jpg'));
%% 5.cal the pixel scale in FFT pattern 换算fft中的像素点尺寸，长宽有区别。要设整张图的pixel大小！
% 每个像素的实际尺寸（纳米）
pixelSize_nm = 0.034;
% 获取图像的尺寸（像素）
imageSize = size(finalManuallyCroppedImage);
% 计算FFT中每个像素的尺寸（1/纳米）
fftPixelSize_per_nm = 1 ./ (pixelSize_nm * imageSize);
disp(['FFT中每个像素的尺寸: ', num2str(fftPixelSize_per_nm), ' 1/nm']);
%% 6.do hanning window fft for croped images 对最终裁剪的图片做fft
% 转换为双精度类型
finalManuallyCroppedImageDouble = im2double(finalManuallyCroppedImage);
% 对裁剪后的图片应用 Hanning 窗口
hanningWindow = hann(size(finalManuallyCroppedImageDouble, 1)) * hann(size(finalManuallyCroppedImageDouble, 2))';
windowedImage = finalManuallyCroppedImageDouble .* hanningWindow;
% 进行二维FFT
fftImage = fft2(windowedImage);
% 将频谱移动至中心
fftImageShifted = fftshift(fftImage);
% 计算频谱的幅度谱
magnitudeSpectrum = abs(fftImageShifted);
% 显示频谱图像
figure;
imshow(log(1 + magnitudeSpectrum), []);
title('Hanning窗口FFT变换后的频谱图');
colormap(hot)
% 归一化
normalizedMagnitudeSpectrum = mat2gray(log(1 + magnitudeSpectrum));
% 应用 colormap
colorMap = hot(256); %改colormap
indexedImage = uint8(ind2rgb(round(normalizedMagnitudeSpectrum * 255) + 1, colorMap) * 255);
% 保存fft到新文件夹
imwrite(indexedImage, fullfile(newFolderPath, 'FFT_of_cropped_image.jpg'));
close
%% 7.plot the 3d image of fft 将fft画成3d图像
% 创建x和y坐标网格
[x_fft3d, y_fft3d] = meshgrid(1:size(magnitudeSpectrum, 2), 1:size(magnitudeSpectrum, 1));
% 创建3D表面图
figure;
surf(x_fft3d, y_fft3d, log(1 + magnitudeSpectrum), 'EdgeColor', 'none');
% 设置图标签
xlabel('X');
ylabel('Y');
zlabel('Log Magnitude Spectrum');
title('FFT image');
% 使用print保存为图像格式
print(fullfile(newFolderPath, '3dFFT_of_cropped_image.png'), '-dpng', '-r300');
close
%% 8.crop the same scale area in vaccum position 在真空中裁剪同样大的区域
% 获取上涨图片裁剪区域的大小
cropWidth = cropPosition(3);
cropHeight = cropPosition(4);
% 新的图片
newImage = imread('G:\++++++data_process\+bigNW_2241_piexl size 0.034nm\output_frames\Aligned 2241 460.0 kx Ceta Camera0000.tif'); % 替换成你的新图片路径
% 使用imrotate函数旋转新的图片
rotatedNewImage = imrotate(newImage, rotationAngle, 'bilinear', 'crop');
% 新的裁剪区域大小
newCropHeight = cropHeight; % 用之前的裁剪高度
newCropWidth = cropWidth;   % 用之前的裁剪宽度
% 显示旋转后的图片
imshow(rotatedNewImage);
title('旋转后的新图片');
% 设置裁剪区域大小的约束函数
positionConstraintFcn = @(position) validateCropSize(position, newCropWidth, newCropHeight);
% 启动交互式裁剪工具，并设置约束函数
h = imrect(gca, 'PositionConstraintFcn', positionConstraintFcn);
% 等待用户手动选择裁剪区域并获取裁剪框的位置信息
wait(h);
newCropPosition = getPosition(h);
% 使用imcrop函数根据用户选择的位置裁剪新的图片
finalManuallyCroppedNewImage = imcrop(rotatedNewImage, newCropPosition);
% 显示最终手动裁剪后的新图片
imshow(finalManuallyCroppedNewImage);
title('最终手动裁剪后的新图片');
% 保存裁剪后的真空图像到新文件夹
imwrite(finalManuallyCroppedNewImage, fullfile(newFolderPath, 'cropped_vaccumImage.jpg'));
close
%% 9.do hanning window fft for cropped vaccum image 对真空区域做同样的fft
% 转换为双精度类型
finalManuallyCroppedNewImageDouble = im2double(finalManuallyCroppedNewImage);
% 对裁剪后的新图片应用 Hanning 窗口
hanningWindowNew = hann(size(finalManuallyCroppedNewImageDouble, 1)) * hann(size(finalManuallyCroppedNewImageDouble, 2))';
windowedImageNew = finalManuallyCroppedNewImageDouble .* hanningWindowNew;
% 进行二维FFT
fftImageNew = fft2(windowedImageNew);
% 将频谱移动至中心
fftImageShiftedNew = fftshift(fftImageNew);
% 计算频谱的幅度谱
magnitudeSpectrumNew = abs(fftImageShiftedNew);
% 显示新裁剪图片的 Hanning 窗口 FFT 变换后的频谱图
figure;
imshow(log(1 + magnitudeSpectrumNew), []);
title('新裁剪图片的 Hanning 窗口 FFT 变换后的频谱图');
colormap(hot)
% 归一化
normalizedMagnitudeSpectrum = mat2gray(log(1 + magnitudeSpectrumNew));
% 应用 colormap
colorMap = hot(256); %改colormap
indexedImage = uint8(ind2rgb(round(normalizedMagnitudeSpectrum * 255) + 1, colorMap) * 255);
% 保存fft到新文件夹
imwrite(indexedImage, fullfile(newFolderPath, 'FFT_of_cropped_vaccumImage.jpg'));
close
%% 10.cal the difference of ffts 计算两个FFT图像的差异，简单看一下
% 计算两个FFT图像的差异
differenceImage = fftImageShifted - fftImageShiftedNew;
% 显示差异图像（使用hot颜色映射）
figure;
% 显示差异图像（使用hot颜色映射），并指定显示范围
imshow(abs(differenceImage), [min(abs(differenceImage(:))), max(abs(differenceImage(:)))]);
title('FFT-FFT of vaccum');
colormap(jet);  % 设置颜色映射为hot
colorbar;       % 显示颜色条
% 归一化差异图像的值范围到 [0, 1]
normalizedDifference = (abs(differenceImage) - min(abs(differenceImage(:)))) / (max(abs(differenceImage(:))) - min(abs(differenceImage(:))));
% 使用jet颜色映射
colormapUsed = hot(256); %修改colormap
colormapImage = uint8(ind2rgb(round(normalizedDifference * 255) + 1, colormapUsed) * 255);
% 保存带有colormap颜色映射的图片
imwrite(colormapImage, fullfile(newFolderPath, 'FFT_difference.jpg'), 'jpg');
close
%% 11.plot the 3d image of fft for vaccum 将真空中的fft画成3d图像
% 创建x和y坐标网格
[x_fft3d, y_fft3d] = meshgrid(1:size(magnitudeSpectrumNew, 2), 1:size(magnitudeSpectrumNew, 1));
% 创建3D表面图
figure;
surf(x_fft3d, y_fft3d, log(1 + magnitudeSpectrumNew), 'EdgeColor', 'none');
% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('Log Magnitude Spectrum');
% 设置图标题
title('FFT image of vaccum');
% 使用print保存为图像格式（如PNG）
print(fullfile(newFolderPath, '3dFFT_of_vaccumImage.png'), '-dpng', '-r300');
close
%% 12.find the brightest piexl in fft of vaccum 找到真空区域fft中的最亮点设为中心点
[maxValue, linearIndex] = max(magnitudeSpectrumNew(:));
[y_center, x_center] = ind2sub(size(magnitudeSpectrumNew), linearIndex);
% 显示新裁剪图片的 Hanning 窗口 FFT 变换后的频谱图
figure;
imshow(log(1 + magnitudeSpectrumNew), []);
hold on;
% 标记真空区域的最亮点
plot(x_center, y_center, 'r*', 'MarkerSize', 10);
title('新裁剪图片的 Hanning 窗口 FFT 变换后的频谱图');
close
%% 13.cal the dis of all piexls to the brightest 计算真空fft中所有像素点到最亮点的距离
% 计算每个像素点到最亮点的距离
distances = sqrt((x_fft3d - x_center).^2 + (y_fft3d - y_center).^2);
% 获取唯一的距离值
uniqueDistances = unique(distances(:));
%% 14.find the brightest piexl per distance 找到所有不同距离下的唯一最亮点
% 找到所有不同距离下的最亮点
brightestPoints = cell(length(uniqueDistances), 1);
% 遍历每个唯一的距离
for i = 1:length(uniqueDistances)
    currentDistance = uniqueDistances(i);
    % 找到当前距离下的最亮点索引
    indicesAtCurrentDistance = find(abs(distances - currentDistance) < 1e-6); % 使用绝对值和一个小的阈值来处理浮点数比较
    % 找到最亮的点索引
    [~, brightestIndex] = max(magnitudeSpectrumNew(indicesAtCurrentDistance));
    brightestIndex = indicesAtCurrentDistance(brightestIndex);
    % 将最亮点的索引转换为坐标
    [y_brightest, x_brightest] = ind2sub(size(distances), brightestIndex);
    % 存储最亮点的坐标
    brightestPoints{i} = [x_brightest, y_brightest];
end
% 移除多余的零元素
brightestPoints = cell2mat(brightestPoints);
% 在图上标记所有不同距离下的最亮点
figure;
imshow(log(1 + magnitudeSpectrumNew), []);
hold on;
plot(brightestPoints(:, 1), brightestPoints(:, 2), 'r*', 'MarkerSize', 10);
title('所有不同距离下的最亮点');
close
%% 15.write the dis to all the same dis piexls 将FFT不同距离的最亮点填写到其他等距离的像素中
% 假设 magnitudeSpectrum 是您要修改的频谱图像
modifiedMagnitudeSpectrum = magnitudeSpectrum;
% 遍历每个唯一的距离
for i = 1:length(uniqueDistances)
    currentDistance = uniqueDistances(i);
    % 找到当前距离下的最亮点索引
    indicesAtCurrentDistance = find(abs(distances - currentDistance) < 1e-6);
    % 找到最亮的点索引
    [~, brightestIndex] = max(magnitudeSpectrumNew(indicesAtCurrentDistance));
    brightestIndex = indicesAtCurrentDistance(brightestIndex);
    % 获取最亮点的值
    brightestValue = magnitudeSpectrumNew(brightestIndex);
    % 将当前距离下的所有像素设置为最亮点的值
    modifiedMagnitudeSpectrum(indicesAtCurrentDistance) = brightestValue;
end
% 显示修改后的频谱图像
figure;
imshow(log(1 + modifiedMagnitudeSpectrum), []);
title('将最亮点的值填充到其他等距离的像素中');
close
%% 16.cal the difference of modified fft again 再次计算两个FFT图像的差异
% 计算两个频谱图像的差异
differenceMagnitudeSpectrum = magnitudeSpectrum - modifiedMagnitudeSpectrum;
% 将所有负值设为零
differenceMagnitudeSpectrum(differenceMagnitudeSpectrum < 0) = 0;
% 显示差异图像
%% 17.save the three fft images 保存三个图
% 创建新的figure，将其子图的大小设置为与屏幕相同
finalFigure = figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(1, 3, 1): 显示原始图像
subplot(1, 3, 1);
imshow(log(1 + magnitudeSpectrum), []);
title('FFT of cropped image');
colorbar;
% subplot(1, 3, 2): 显示修改后的图像
subplot(1, 3, 2);
imshow(log(1 + modifiedMagnitudeSpectrum), []);
title('modified FFT of vaccum image');
colorbar;
% subplot(1, 3, 3): 显示两者的逐像素差异
subplot(1, 3, 3);
imshow(log(1 + differenceMagnitudeSpectrum), []);
title('modified FFT difference');
colorbar;
% 设置颜色映射
colormap(parula);
% 将figure放到最大分辨率
% 构建要保存的文件的完整路径
filename = 'final_fft.jpg';
fullFilePath = fullfile(newFolderPath, filename);
% 保存图像为 JPG 文件，直接使用 print 函数
print(finalFigure, fullFilePath, '-djpeg', '-r400'); % '-djpeg' 指定保存为 JPEG 格式，'-r300' 指定分辨率为 300 dpi
close
%% 18.instead of spectrum 参数换回去，计算距离与夹角
magnitudeSpectrum = differenceMagnitudeSpectrum;
%% 19.cal the center in fft image 计算fft图的中心pixel
% 获取频谱图的尺寸
[m, n] = size(magnitudeSpectrum);
% 计算中心点的坐标
centerX = ceil(n / 2);
centerY = ceil(m / 2);
disp(['FFT中心点坐标：(', num2str(centerX), ', ', num2str(centerY), ')']);
%% ***找fft的实际中心，如果fft是偶数需要调整中心坐标(如果下一步不行就用这一步)
% 获取频谱图的尺寸
[m, n] = size(magnitudeSpectrum);

% 计算中心点的坐标
centerX = ceil(n / 2);
centerY = ceil(m / 2);

% 如果频谱图的尺寸是偶数，则需要调整中心坐标
if mod(n, 2) == 0
    centerX = centerX + 1;
end

if mod(m, 2) == 0
    centerY = centerY + 1;
end
disp(['FFT中心点坐标：(', num2str(centerX), ', ', num2str(centerY), ')']);
%% 20. cal the position of the most bright point 计算fft中最亮的pixel，并设为原点
% 找到频谱图中的最大值及其坐标
[maxValue, maxIndex] = max(magnitudeSpectrum(:));
% 将线性索引(maxIndex)转换为行列索引
[row, col] = ind2sub(size(magnitudeSpectrum), maxIndex);
disp(['最亮点的幅度值：', num2str(maxValue)]);
disp(['最亮点的坐标：(', num2str(col), ', ', num2str(row), ')']);
%让原点等于最亮的点位置
centerX=col;
centerY=row;
%% 21.delete the piexls by lattice limit 通过最小晶面间距，扣去中心点附近的pixels，要设参数！
% 假设每个像素的实际大小
pixelSizeInMillimetersVertical = fftPixelSize_per_nm(:,1);
pixelSizeInMillimetersHorizontal = fftPixelSize_per_nm(:,2);
% 计算原点附近的物理范围
centerRadiusInMillimeters = 2.479/2;  % 1/nm除以2，设置为Bi最小晶面间距倒数的大小,这里用的是003晶面2.479 1/nm
% 将物理尺度转换为像素尺度
centerRadiusVertical = round(centerRadiusInMillimeters / pixelSizeInMillimetersVertical);
centerRadiusHorizontal = round(centerRadiusInMillimeters / pixelSizeInMillimetersHorizontal);
% 将原点附近的点的幅度设置为零
magnitudeSpectrum(centerY - centerRadiusVertical:centerY + centerRadiusVertical, ...
                   centerX - centerRadiusHorizontal:centerX + centerRadiusHorizontal) = 0;
% 显示更新后的频谱图像
figure;
imshow(log(1 + magnitudeSpectrum), []);
title('去除原点附近的点后的频谱图');
close
%% 不用运行 delete the pixel by the area around the center 和上面一样，用范围扣掉原点附近的一定范围内的点，要设范围！
% 假设每个像素的实际大小（垂直方向，水平方向）
pixelSizeInMillimetersVertical = fftPixelSize_per_nm(:,1);
pixelSizeInMillimetersHorizontal = fftPixelSize_per_nm(:,2);
% 设置原点附近的点的半径，可以根据需要调整
centerRadius = 5;
% 将原点附近的点的幅度设置为零
magnitudeSpectrum(centerY - centerRadius:centerY + centerRadius, centerX - centerRadius:centerX + centerRadius) = 0;
% 显示更新后的频谱图像
figure;
imshow(log(1 + magnitudeSpectrum), []);
title('去除原点附近的点后的频谱图');
close
%% 22.find other most bright pixels in fft 找其他比较亮的点。要设阈值！
% 寻找频谱图中的峰值
threshold = 0.3;  % 设置阈值，可以根据需要调整
peaks = imregionalmax(magnitudeSpectrum, 4) & (magnitudeSpectrum > threshold * max(magnitudeSpectrum(:)));
% 获取峰值的坐标
[peakRows, peakCols] = find(peaks);
% 显示处理后的频谱图和标记峰值
figure;
imshow(log(1 + magnitudeSpectrum), []);
hold on;
plot(peakCols, peakRows, 'r*');
title('处理后的频谱图中的亮点');
%close
%% 23.mark the distance and angle of other pixels 计算距离与垂直线之间的夹角。注意pixel大小！
% 假设每个像素的实际大小（垂直方向，水平方向）
pixelSizeInMillimetersVertical = fftPixelSize_per_nm(:,1);
pixelSizeInMillimetersHorizontal = fftPixelSize_per_nm(:,2);
% 计算每个点与原点之间的距离（以像素为单位）
distances = sqrt(((peakRows - centerY) * pixelSizeInMillimetersVertical).^2 + ((peakCols - centerX) * pixelSizeInMillimetersHorizontal).^2);
% 计算每个点与原点之间的夹角（弧度）
anglesRad = atan2((peakCols - centerX) * pixelSizeInMillimetersHorizontal, (peakRows - centerY) * pixelSizeInMillimetersVertical);
% 将弧度转换为角
anglesDegrees = rad2deg(anglesRad);
% 显示处理后的频谱图和标记峰值、距离、夹角、连线和垂直线
figure;
imshow(log(1 + magnitudeSpectrum), []);
hold on;
plot(peakCols, peakRows, 'r*');
% 在图中标记亮点与原点的距离、夹角和连线
for i = 1:numel(peakRows)
    % 画出每个点与原点的连线
    plot([centerX, peakCols(i)], [centerY, peakRows(i)], 'g-');
    % 画出垂直线
    plot([centerX, centerX], [centerY, peakRows(i)], 'w-');
    % 在图中标记距离和夹角，保留小数点后两位
    text(peakCols(i), peakRows(i), [sprintf('%.1f', distances(i)), '1/nm, ', sprintf('%.1f', anglesDegrees(i)), '°'], 'Color', 'yellow');
end
title('处理后的频谱图中的亮点、距离、夹角、连线和垂直线');
%% 24.save the fft and dis. 保留fft和量距离的fig文件
figure;
% 显示处理后的频谱图
imshow(log(1 + magnitudeSpectrum), []);
hold on;
% 标记峰值位置
plot(peakCols, peakRows, 'r*');
% 绘制线条、垂直线，并添加只显示距离的文本注释
textSeparation = -2; % 根据需要调整这个值
for i = 1:numel(peakRows)
    plot([centerX, peakCols(i)], [centerY, peakRows(i)], 'g-');
    plot([centerX, centerX], [centerY, peakRows(i)], 'w-');
    % 计算文本位置
    textX = peakCols(i) + textSeparation;
    textY = peakRows(i) + textSeparation;
    % 显示距离值
    text(textX, textY, [sprintf('%.1f', distances(i)), '1/nm'], 'Color', 'yellow');
end
% 添加标题
title('处理后的频谱图中的亮点、距离、连线和垂直线');
% 保存图形为.fig格式（高分辨率）
savefig(fullfile(newFolderPath, 'modified_fft_w_dis.fig'));
close
%% 25.save fft and angles 保留fft和量的角度的fig文件
% 创建一个图形
figure;
% 显示处理后的频谱图
imshow(log(1 + magnitudeSpectrum), []);
hold on;
% 标记峰值位置
plot(peakCols, peakRows, 'r*');
% 绘制线条、垂直线，并添加只显示角度的文本注释
textSeparation = -2; % 根据需要调整这个值
for i = 1:numel(peakRows)
    plot([centerX, peakCols(i)], [centerY, peakRows(i)], 'g-');
    plot([centerX, centerX], [centerY, peakRows(i)], 'w-');
    % 计算文本位置
    textX = peakCols(i) + textSeparation;
    textY = peakRows(i) + textSeparation;
    % 显示角度值
    text(textX, textY, [sprintf('%.1f', anglesDegrees(i)), '°'], 'Color', 'yellow');
end
% 添加标题
title('处理后的频谱图中的亮点、夹角、连线和垂直线');
% 保存图形为.fig格式（高分辨率）
savefig(fullfile(newFolderPath, 'modified_fft_w_angles.fig'));
close
%% 26.acquire the row, distance and angle 将行数，距离，角度写进output中
% 获取完整的文件名（包括路径）
fullFileName = imagePath
% 使用fileparts获取文件名和扩展名
[~, fileName, ~] = fileparts(fullFileName);
% 获取文件名的后四位
fileSuffix = fileName(end-3:end);
% 计算每个点与原点之间的距离（以毫米为单位）
distances = sqrt(((peakCols - centerX) * pixelSizeInMillimetersHorizontal).^2 + ((peakRows - centerY) * pixelSizeInMillimetersVertical).^2);
% 取前一半的数据
halfIndex = ceil(numel(distances)/2);
distances = distances(1:halfIndex);
anglesDegrees = anglesDegrees(1:halfIndex);
% 排序距离，并根据排序结果重新排列角度
[sortedDistances, sortedIndices] = sort(distances);
sortedAnglesDegrees = anglesDegrees(sortedIndices);
% 对距离和角度保留两位小数
sortedDistances = round(sortedDistances, 2);
sortedAnglesDegrees = round(sortedAnglesDegrees, 1);
% 将舍入后的距离和角度分配给表格
output_dis_angle(str2double(fileSuffix), :) = {sortedDistances', sortedAnglesDegrees'};
% 创建一个表格，包含 Distance 和 Angle 列
output_dis_angle(str2double(fileSuffix), :) = {sortedDistances', sortedAnglesDegrees'};
% 在 MATLAB 工作区中显示表格
disp(output_dis_angle);
%% 27.保存output为mat文件，下次可以继续将参数写进来
outputFilePath = 'G:\++++++data_process\+bigNW_2241_piexl size 0.034nm\output_frames/output_dis_angle.mat';
save(outputFilePath, 'output_dis_angle');
%%

