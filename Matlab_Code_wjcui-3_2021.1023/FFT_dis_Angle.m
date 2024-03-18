%% 0.clear
clc
clear
load('G:\++++++data_process\+bigNW_2241_piexl size 0.034nm\output_frames\output_dis_angle.mat')
%% 1.read image with tif format ����tifͼƬ
imagePath = 'G:\++++++data_process\+bigNW_2241_piexl size 0.034nm\output_frames\Aligned 2241 460.0 kx Ceta Camera0553.tif'; % ���滻 'your_image.tif' Ϊʵ��ͼƬ�ļ���·��
imageData = imread(imagePath); % ʹ��imread������ȡtifͼƬ
%% 2.crop the import image ���²ü�һ��ͼƬ
imageData = imread(imagePath); % ʹ��imread������ȡtifͼƬ
imshow(imageData);
title('�����ͼƬ'); % ��ʾ�����ͼƬ
croppedImage = imcrop; % ʹ��imcrop�����ֶ��ü�ͼ��
figure;  
imshow(croppedImage);
title('�ü����ͼƬ'); % ��ʾ�ü����ͼ��
close
%% 3.rotate the image ��תͼƬ������תˮƽ���ٴβü���Ҫ����ת�Ƕȣ�
% ָ��ͼƬ��ת�ĽǶȲ���
rotationAngle = 74;
% ʹ��imrotate������תͼƬ
rotatedImage = imrotate(croppedImage, rotationAngle, 'bilinear', 'crop');
% ��ʾ��ת���ͼƬ
imshow(rotatedImage);
title('��ת���ͼƬ');
% ��������ʽ�ü�����
h = imrect;
% �ȴ��û��ֶ�ѡ��ü����򲢻�ȡ�ü����λ����Ϣ
wait(h);
cropPosition = getPosition(h);
% ʹ��imcrop���������û�ѡ���λ�òü���ת���ͼƬ
finalManuallyCroppedImage = imcrop(rotatedImage, cropPosition);
% ��ʾ�����ֶ��ü����ͼƬ
imshow(finalManuallyCroppedImage);
title('�����ֶ��ü����ͼƬ');
close
%% 4.save the final corpped image �������ղü���ͼƬ
% ��ȡ��һ��Ŀ¼·��
parentFolder = fileparts(imagePath);
% ��ȡԭʼ�ļ����ĺ���λ
[~, filename, ~] = fileparts(imagePath);
suffix = filename(end-3:end);
% �������ļ��У��ļ�����Ϊԭʼ�ļ����ĺ���λ
newFolderName = suffix;
newFolderPath = fullfile(parentFolder, newFolderName);
% ����ļ��в����ڣ��򴴽�
if ~exist(newFolderPath, 'dir')
    mkdir(newFolderPath);
end
% ����ü����ͼ�����ļ���
imwrite(finalManuallyCroppedImage, fullfile(newFolderPath, 'cropped_image.jpg'));
%% 5.cal the pixel scale in FFT pattern ����fft�е����ص�ߴ磬����������Ҫ������ͼ��pixel��С��
% ÿ�����ص�ʵ�ʳߴ磨���ף�
pixelSize_nm = 0.034;
% ��ȡͼ��ĳߴ磨���أ�
imageSize = size(finalManuallyCroppedImage);
% ����FFT��ÿ�����صĳߴ磨1/���ף�
fftPixelSize_per_nm = 1 ./ (pixelSize_nm * imageSize);
disp(['FFT��ÿ�����صĳߴ�: ', num2str(fftPixelSize_per_nm), ' 1/nm']);
%% 6.do hanning window fft for croped images �����ղü���ͼƬ��fft
% ת��Ϊ˫��������
finalManuallyCroppedImageDouble = im2double(finalManuallyCroppedImage);
% �Բü����ͼƬӦ�� Hanning ����
hanningWindow = hann(size(finalManuallyCroppedImageDouble, 1)) * hann(size(finalManuallyCroppedImageDouble, 2))';
windowedImage = finalManuallyCroppedImageDouble .* hanningWindow;
% ���ж�άFFT
fftImage = fft2(windowedImage);
% ��Ƶ���ƶ�������
fftImageShifted = fftshift(fftImage);
% ����Ƶ�׵ķ�����
magnitudeSpectrum = abs(fftImageShifted);
% ��ʾƵ��ͼ��
figure;
imshow(log(1 + magnitudeSpectrum), []);
title('Hanning����FFT�任���Ƶ��ͼ');
colormap(hot)
% ��һ��
normalizedMagnitudeSpectrum = mat2gray(log(1 + magnitudeSpectrum));
% Ӧ�� colormap
colorMap = hot(256); %��colormap
indexedImage = uint8(ind2rgb(round(normalizedMagnitudeSpectrum * 255) + 1, colorMap) * 255);
% ����fft�����ļ���
imwrite(indexedImage, fullfile(newFolderPath, 'FFT_of_cropped_image.jpg'));
close
%% 7.plot the 3d image of fft ��fft����3dͼ��
% ����x��y��������
[x_fft3d, y_fft3d] = meshgrid(1:size(magnitudeSpectrum, 2), 1:size(magnitudeSpectrum, 1));
% ����3D����ͼ
figure;
surf(x_fft3d, y_fft3d, log(1 + magnitudeSpectrum), 'EdgeColor', 'none');
% ����ͼ��ǩ
xlabel('X');
ylabel('Y');
zlabel('Log Magnitude Spectrum');
title('FFT image');
% ʹ��print����Ϊͼ���ʽ
print(fullfile(newFolderPath, '3dFFT_of_cropped_image.png'), '-dpng', '-r300');
close
%% 8.crop the same scale area in vaccum position ������вü�ͬ���������
% ��ȡ����ͼƬ�ü�����Ĵ�С
cropWidth = cropPosition(3);
cropHeight = cropPosition(4);
% �µ�ͼƬ
newImage = imread('G:\++++++data_process\+bigNW_2241_piexl size 0.034nm\output_frames\Aligned 2241 460.0 kx Ceta Camera0000.tif'); % �滻�������ͼƬ·��
% ʹ��imrotate������ת�µ�ͼƬ
rotatedNewImage = imrotate(newImage, rotationAngle, 'bilinear', 'crop');
% �µĲü������С
newCropHeight = cropHeight; % ��֮ǰ�Ĳü��߶�
newCropWidth = cropWidth;   % ��֮ǰ�Ĳü����
% ��ʾ��ת���ͼƬ
imshow(rotatedNewImage);
title('��ת�����ͼƬ');
% ���òü������С��Լ������
positionConstraintFcn = @(position) validateCropSize(position, newCropWidth, newCropHeight);
% ��������ʽ�ü����ߣ�������Լ������
h = imrect(gca, 'PositionConstraintFcn', positionConstraintFcn);
% �ȴ��û��ֶ�ѡ��ü����򲢻�ȡ�ü����λ����Ϣ
wait(h);
newCropPosition = getPosition(h);
% ʹ��imcrop���������û�ѡ���λ�òü��µ�ͼƬ
finalManuallyCroppedNewImage = imcrop(rotatedNewImage, newCropPosition);
% ��ʾ�����ֶ��ü������ͼƬ
imshow(finalManuallyCroppedNewImage);
title('�����ֶ��ü������ͼƬ');
% ����ü�������ͼ�����ļ���
imwrite(finalManuallyCroppedNewImage, fullfile(newFolderPath, 'cropped_vaccumImage.jpg'));
close
%% 9.do hanning window fft for cropped vaccum image �����������ͬ����fft
% ת��Ϊ˫��������
finalManuallyCroppedNewImageDouble = im2double(finalManuallyCroppedNewImage);
% �Բü������ͼƬӦ�� Hanning ����
hanningWindowNew = hann(size(finalManuallyCroppedNewImageDouble, 1)) * hann(size(finalManuallyCroppedNewImageDouble, 2))';
windowedImageNew = finalManuallyCroppedNewImageDouble .* hanningWindowNew;
% ���ж�άFFT
fftImageNew = fft2(windowedImageNew);
% ��Ƶ���ƶ�������
fftImageShiftedNew = fftshift(fftImageNew);
% ����Ƶ�׵ķ�����
magnitudeSpectrumNew = abs(fftImageShiftedNew);
% ��ʾ�²ü�ͼƬ�� Hanning ���� FFT �任���Ƶ��ͼ
figure;
imshow(log(1 + magnitudeSpectrumNew), []);
title('�²ü�ͼƬ�� Hanning ���� FFT �任���Ƶ��ͼ');
colormap(hot)
% ��һ��
normalizedMagnitudeSpectrum = mat2gray(log(1 + magnitudeSpectrumNew));
% Ӧ�� colormap
colorMap = hot(256); %��colormap
indexedImage = uint8(ind2rgb(round(normalizedMagnitudeSpectrum * 255) + 1, colorMap) * 255);
% ����fft�����ļ���
imwrite(indexedImage, fullfile(newFolderPath, 'FFT_of_cropped_vaccumImage.jpg'));
close
%% 10.cal the difference of ffts ��������FFTͼ��Ĳ��죬�򵥿�һ��
% ��������FFTͼ��Ĳ���
differenceImage = fftImageShifted - fftImageShiftedNew;
% ��ʾ����ͼ��ʹ��hot��ɫӳ�䣩
figure;
% ��ʾ����ͼ��ʹ��hot��ɫӳ�䣩����ָ����ʾ��Χ
imshow(abs(differenceImage), [min(abs(differenceImage(:))), max(abs(differenceImage(:)))]);
title('FFT-FFT of vaccum');
colormap(jet);  % ������ɫӳ��Ϊhot
colorbar;       % ��ʾ��ɫ��
% ��һ������ͼ���ֵ��Χ�� [0, 1]
normalizedDifference = (abs(differenceImage) - min(abs(differenceImage(:)))) / (max(abs(differenceImage(:))) - min(abs(differenceImage(:))));
% ʹ��jet��ɫӳ��
colormapUsed = hot(256); %�޸�colormap
colormapImage = uint8(ind2rgb(round(normalizedDifference * 255) + 1, colormapUsed) * 255);
% �������colormap��ɫӳ���ͼƬ
imwrite(colormapImage, fullfile(newFolderPath, 'FFT_difference.jpg'), 'jpg');
close
%% 11.plot the 3d image of fft for vaccum ������е�fft����3dͼ��
% ����x��y��������
[x_fft3d, y_fft3d] = meshgrid(1:size(magnitudeSpectrumNew, 2), 1:size(magnitudeSpectrumNew, 1));
% ����3D����ͼ
figure;
surf(x_fft3d, y_fft3d, log(1 + magnitudeSpectrumNew), 'EdgeColor', 'none');
% �����������ǩ
xlabel('X');
ylabel('Y');
zlabel('Log Magnitude Spectrum');
% ����ͼ����
title('FFT image of vaccum');
% ʹ��print����Ϊͼ���ʽ����PNG��
print(fullfile(newFolderPath, '3dFFT_of_vaccumImage.png'), '-dpng', '-r300');
close
%% 12.find the brightest piexl in fft of vaccum �ҵ��������fft�е���������Ϊ���ĵ�
[maxValue, linearIndex] = max(magnitudeSpectrumNew(:));
[y_center, x_center] = ind2sub(size(magnitudeSpectrumNew), linearIndex);
% ��ʾ�²ü�ͼƬ�� Hanning ���� FFT �任���Ƶ��ͼ
figure;
imshow(log(1 + magnitudeSpectrumNew), []);
hold on;
% �����������������
plot(x_center, y_center, 'r*', 'MarkerSize', 10);
title('�²ü�ͼƬ�� Hanning ���� FFT �任���Ƶ��ͼ');
close
%% 13.cal the dis of all piexls to the brightest �������fft���������ص㵽������ľ���
% ����ÿ�����ص㵽������ľ���
distances = sqrt((x_fft3d - x_center).^2 + (y_fft3d - y_center).^2);
% ��ȡΨһ�ľ���ֵ
uniqueDistances = unique(distances(:));
%% 14.find the brightest piexl per distance �ҵ����в�ͬ�����µ�Ψһ������
% �ҵ����в�ͬ�����µ�������
brightestPoints = cell(length(uniqueDistances), 1);
% ����ÿ��Ψһ�ľ���
for i = 1:length(uniqueDistances)
    currentDistance = uniqueDistances(i);
    % �ҵ���ǰ�����µ�����������
    indicesAtCurrentDistance = find(abs(distances - currentDistance) < 1e-6); % ʹ�þ���ֵ��һ��С����ֵ�����������Ƚ�
    % �ҵ������ĵ�����
    [~, brightestIndex] = max(magnitudeSpectrumNew(indicesAtCurrentDistance));
    brightestIndex = indicesAtCurrentDistance(brightestIndex);
    % �������������ת��Ϊ����
    [y_brightest, x_brightest] = ind2sub(size(distances), brightestIndex);
    % �洢�����������
    brightestPoints{i} = [x_brightest, y_brightest];
end
% �Ƴ��������Ԫ��
brightestPoints = cell2mat(brightestPoints);
% ��ͼ�ϱ�����в�ͬ�����µ�������
figure;
imshow(log(1 + magnitudeSpectrumNew), []);
hold on;
plot(brightestPoints(:, 1), brightestPoints(:, 2), 'r*', 'MarkerSize', 10);
title('���в�ͬ�����µ�������');
close
%% 15.write the dis to all the same dis piexls ��FFT��ͬ�������������д�������Ⱦ����������
% ���� magnitudeSpectrum ����Ҫ�޸ĵ�Ƶ��ͼ��
modifiedMagnitudeSpectrum = magnitudeSpectrum;
% ����ÿ��Ψһ�ľ���
for i = 1:length(uniqueDistances)
    currentDistance = uniqueDistances(i);
    % �ҵ���ǰ�����µ�����������
    indicesAtCurrentDistance = find(abs(distances - currentDistance) < 1e-6);
    % �ҵ������ĵ�����
    [~, brightestIndex] = max(magnitudeSpectrumNew(indicesAtCurrentDistance));
    brightestIndex = indicesAtCurrentDistance(brightestIndex);
    % ��ȡ�������ֵ
    brightestValue = magnitudeSpectrumNew(brightestIndex);
    % ����ǰ�����µ�������������Ϊ�������ֵ
    modifiedMagnitudeSpectrum(indicesAtCurrentDistance) = brightestValue;
end
% ��ʾ�޸ĺ��Ƶ��ͼ��
figure;
imshow(log(1 + modifiedMagnitudeSpectrum), []);
title('���������ֵ��䵽�����Ⱦ����������');
close
%% 16.cal the difference of modified fft again �ٴμ�������FFTͼ��Ĳ���
% ��������Ƶ��ͼ��Ĳ���
differenceMagnitudeSpectrum = magnitudeSpectrum - modifiedMagnitudeSpectrum;
% �����и�ֵ��Ϊ��
differenceMagnitudeSpectrum(differenceMagnitudeSpectrum < 0) = 0;
% ��ʾ����ͼ��
%% 17.save the three fft images ��������ͼ
% �����µ�figure��������ͼ�Ĵ�С����Ϊ����Ļ��ͬ
finalFigure = figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(1, 3, 1): ��ʾԭʼͼ��
subplot(1, 3, 1);
imshow(log(1 + magnitudeSpectrum), []);
title('FFT of cropped image');
colorbar;
% subplot(1, 3, 2): ��ʾ�޸ĺ��ͼ��
subplot(1, 3, 2);
imshow(log(1 + modifiedMagnitudeSpectrum), []);
title('modified FFT of vaccum image');
colorbar;
% subplot(1, 3, 3): ��ʾ���ߵ������ز���
subplot(1, 3, 3);
imshow(log(1 + differenceMagnitudeSpectrum), []);
title('modified FFT difference');
colorbar;
% ������ɫӳ��
colormap(parula);
% ��figure�ŵ����ֱ���
% ����Ҫ������ļ�������·��
filename = 'final_fft.jpg';
fullFilePath = fullfile(newFolderPath, filename);
% ����ͼ��Ϊ JPG �ļ���ֱ��ʹ�� print ����
print(finalFigure, fullFilePath, '-djpeg', '-r400'); % '-djpeg' ָ������Ϊ JPEG ��ʽ��'-r300' ָ���ֱ���Ϊ 300 dpi
close
%% 18.instead of spectrum ��������ȥ�����������н�
magnitudeSpectrum = differenceMagnitudeSpectrum;
%% 19.cal the center in fft image ����fftͼ������pixel
% ��ȡƵ��ͼ�ĳߴ�
[m, n] = size(magnitudeSpectrum);
% �������ĵ������
centerX = ceil(n / 2);
centerY = ceil(m / 2);
disp(['FFT���ĵ����꣺(', num2str(centerX), ', ', num2str(centerY), ')']);
%% ***��fft��ʵ�����ģ����fft��ż����Ҫ������������(�����һ�����о�����һ��)
% ��ȡƵ��ͼ�ĳߴ�
[m, n] = size(magnitudeSpectrum);

% �������ĵ������
centerX = ceil(n / 2);
centerY = ceil(m / 2);

% ���Ƶ��ͼ�ĳߴ���ż��������Ҫ������������
if mod(n, 2) == 0
    centerX = centerX + 1;
end

if mod(m, 2) == 0
    centerY = centerY + 1;
end
disp(['FFT���ĵ����꣺(', num2str(centerX), ', ', num2str(centerY), ')']);
%% 20. cal the position of the most bright point ����fft��������pixel������Ϊԭ��
% �ҵ�Ƶ��ͼ�е����ֵ��������
[maxValue, maxIndex] = max(magnitudeSpectrum(:));
% ����������(maxIndex)ת��Ϊ��������
[row, col] = ind2sub(size(magnitudeSpectrum), maxIndex);
disp(['������ķ���ֵ��', num2str(maxValue)]);
disp(['����������꣺(', num2str(col), ', ', num2str(row), ')']);
%��ԭ����������ĵ�λ��
centerX=col;
centerY=row;
%% 21.delete the piexls by lattice limit ͨ����С�����࣬��ȥ���ĵ㸽����pixels��Ҫ�������
% ����ÿ�����ص�ʵ�ʴ�С
pixelSizeInMillimetersVertical = fftPixelSize_per_nm(:,1);
pixelSizeInMillimetersHorizontal = fftPixelSize_per_nm(:,2);
% ����ԭ�㸽��������Χ
centerRadiusInMillimeters = 2.479/2;  % 1/nm����2������ΪBi��С�����൹���Ĵ�С,�����õ���003����2.479 1/nm
% ������߶�ת��Ϊ���س߶�
centerRadiusVertical = round(centerRadiusInMillimeters / pixelSizeInMillimetersVertical);
centerRadiusHorizontal = round(centerRadiusInMillimeters / pixelSizeInMillimetersHorizontal);
% ��ԭ�㸽���ĵ�ķ�������Ϊ��
magnitudeSpectrum(centerY - centerRadiusVertical:centerY + centerRadiusVertical, ...
                   centerX - centerRadiusHorizontal:centerX + centerRadiusHorizontal) = 0;
% ��ʾ���º��Ƶ��ͼ��
figure;
imshow(log(1 + magnitudeSpectrum), []);
title('ȥ��ԭ�㸽���ĵ���Ƶ��ͼ');
close
%% �������� delete the pixel by the area around the center ������һ�����÷�Χ�۵�ԭ�㸽����һ����Χ�ڵĵ㣬Ҫ�跶Χ��
% ����ÿ�����ص�ʵ�ʴ�С����ֱ����ˮƽ����
pixelSizeInMillimetersVertical = fftPixelSize_per_nm(:,1);
pixelSizeInMillimetersHorizontal = fftPixelSize_per_nm(:,2);
% ����ԭ�㸽���ĵ�İ뾶�����Ը�����Ҫ����
centerRadius = 5;
% ��ԭ�㸽���ĵ�ķ�������Ϊ��
magnitudeSpectrum(centerY - centerRadius:centerY + centerRadius, centerX - centerRadius:centerX + centerRadius) = 0;
% ��ʾ���º��Ƶ��ͼ��
figure;
imshow(log(1 + magnitudeSpectrum), []);
title('ȥ��ԭ�㸽���ĵ���Ƶ��ͼ');
close
%% 22.find other most bright pixels in fft �������Ƚ����ĵ㡣Ҫ����ֵ��
% Ѱ��Ƶ��ͼ�еķ�ֵ
threshold = 0.3;  % ������ֵ�����Ը�����Ҫ����
peaks = imregionalmax(magnitudeSpectrum, 4) & (magnitudeSpectrum > threshold * max(magnitudeSpectrum(:)));
% ��ȡ��ֵ������
[peakRows, peakCols] = find(peaks);
% ��ʾ������Ƶ��ͼ�ͱ�Ƿ�ֵ
figure;
imshow(log(1 + magnitudeSpectrum), []);
hold on;
plot(peakCols, peakRows, 'r*');
title('������Ƶ��ͼ�е�����');
%close
%% 23.mark the distance and angle of other pixels ��������봹ֱ��֮��ļнǡ�ע��pixel��С��
% ����ÿ�����ص�ʵ�ʴ�С����ֱ����ˮƽ����
pixelSizeInMillimetersVertical = fftPixelSize_per_nm(:,1);
pixelSizeInMillimetersHorizontal = fftPixelSize_per_nm(:,2);
% ����ÿ������ԭ��֮��ľ��루������Ϊ��λ��
distances = sqrt(((peakRows - centerY) * pixelSizeInMillimetersVertical).^2 + ((peakCols - centerX) * pixelSizeInMillimetersHorizontal).^2);
% ����ÿ������ԭ��֮��ļнǣ����ȣ�
anglesRad = atan2((peakCols - centerX) * pixelSizeInMillimetersHorizontal, (peakRows - centerY) * pixelSizeInMillimetersVertical);
% ������ת��Ϊ��
anglesDegrees = rad2deg(anglesRad);
% ��ʾ������Ƶ��ͼ�ͱ�Ƿ�ֵ�����롢�нǡ����ߺʹ�ֱ��
figure;
imshow(log(1 + magnitudeSpectrum), []);
hold on;
plot(peakCols, peakRows, 'r*');
% ��ͼ�б��������ԭ��ľ��롢�нǺ�����
for i = 1:numel(peakRows)
    % ����ÿ������ԭ�������
    plot([centerX, peakCols(i)], [centerY, peakRows(i)], 'g-');
    % ������ֱ��
    plot([centerX, centerX], [centerY, peakRows(i)], 'w-');
    % ��ͼ�б�Ǿ���ͼнǣ�����С�������λ
    text(peakCols(i), peakRows(i), [sprintf('%.1f', distances(i)), '1/nm, ', sprintf('%.1f', anglesDegrees(i)), '��'], 'Color', 'yellow');
end
title('������Ƶ��ͼ�е����㡢���롢�нǡ����ߺʹ�ֱ��');
%% 24.save the fft and dis. ����fft���������fig�ļ�
figure;
% ��ʾ������Ƶ��ͼ
imshow(log(1 + magnitudeSpectrum), []);
hold on;
% ��Ƿ�ֵλ��
plot(peakCols, peakRows, 'r*');
% ������������ֱ�ߣ������ֻ��ʾ������ı�ע��
textSeparation = -2; % ������Ҫ�������ֵ
for i = 1:numel(peakRows)
    plot([centerX, peakCols(i)], [centerY, peakRows(i)], 'g-');
    plot([centerX, centerX], [centerY, peakRows(i)], 'w-');
    % �����ı�λ��
    textX = peakCols(i) + textSeparation;
    textY = peakRows(i) + textSeparation;
    % ��ʾ����ֵ
    text(textX, textY, [sprintf('%.1f', distances(i)), '1/nm'], 'Color', 'yellow');
end
% ��ӱ���
title('������Ƶ��ͼ�е����㡢���롢���ߺʹ�ֱ��');
% ����ͼ��Ϊ.fig��ʽ���߷ֱ��ʣ�
savefig(fullfile(newFolderPath, 'modified_fft_w_dis.fig'));
close
%% 25.save fft and angles ����fft�����ĽǶȵ�fig�ļ�
% ����һ��ͼ��
figure;
% ��ʾ������Ƶ��ͼ
imshow(log(1 + magnitudeSpectrum), []);
hold on;
% ��Ƿ�ֵλ��
plot(peakCols, peakRows, 'r*');
% ������������ֱ�ߣ������ֻ��ʾ�Ƕȵ��ı�ע��
textSeparation = -2; % ������Ҫ�������ֵ
for i = 1:numel(peakRows)
    plot([centerX, peakCols(i)], [centerY, peakRows(i)], 'g-');
    plot([centerX, centerX], [centerY, peakRows(i)], 'w-');
    % �����ı�λ��
    textX = peakCols(i) + textSeparation;
    textY = peakRows(i) + textSeparation;
    % ��ʾ�Ƕ�ֵ
    text(textX, textY, [sprintf('%.1f', anglesDegrees(i)), '��'], 'Color', 'yellow');
end
% ��ӱ���
title('������Ƶ��ͼ�е����㡢�нǡ����ߺʹ�ֱ��');
% ����ͼ��Ϊ.fig��ʽ���߷ֱ��ʣ�
savefig(fullfile(newFolderPath, 'modified_fft_w_angles.fig'));
close
%% 26.acquire the row, distance and angle �����������룬�Ƕ�д��output��
% ��ȡ�������ļ���������·����
fullFileName = imagePath
% ʹ��fileparts��ȡ�ļ�������չ��
[~, fileName, ~] = fileparts(fullFileName);
% ��ȡ�ļ����ĺ���λ
fileSuffix = fileName(end-3:end);
% ����ÿ������ԭ��֮��ľ��루�Ժ���Ϊ��λ��
distances = sqrt(((peakCols - centerX) * pixelSizeInMillimetersHorizontal).^2 + ((peakRows - centerY) * pixelSizeInMillimetersVertical).^2);
% ȡǰһ�������
halfIndex = ceil(numel(distances)/2);
distances = distances(1:halfIndex);
anglesDegrees = anglesDegrees(1:halfIndex);
% ������룬�������������������нǶ�
[sortedDistances, sortedIndices] = sort(distances);
sortedAnglesDegrees = anglesDegrees(sortedIndices);
% �Ծ���ͽǶȱ�����λС��
sortedDistances = round(sortedDistances, 2);
sortedAnglesDegrees = round(sortedAnglesDegrees, 1);
% �������ľ���ͽǶȷ�������
output_dis_angle(str2double(fileSuffix), :) = {sortedDistances', sortedAnglesDegrees'};
% ����һ����񣬰��� Distance �� Angle ��
output_dis_angle(str2double(fileSuffix), :) = {sortedDistances', sortedAnglesDegrees'};
% �� MATLAB ����������ʾ���
disp(output_dis_angle);
%% 27.����outputΪmat�ļ����´ο��Լ���������д����
outputFilePath = 'G:\++++++data_process\+bigNW_2241_piexl size 0.034nm\output_frames/output_dis_angle.mat';
save(outputFilePath, 'output_dis_angle');
%%

