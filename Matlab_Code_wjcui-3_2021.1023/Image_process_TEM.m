%% 0.clear
clc
clear
%% 1.read image with tif format
imagePath = 'D:\++++++Data\Matlab Data\Exp_Data\Bi liquid\+relat fft and image\791.tif'; % ���滻 'your_image.tif' Ϊʵ��ͼƬ�ļ���·��
imageData = imread(imagePath); % ʹ��imread������ȡtifͼƬ
croppedImageData = imageData;
%% ++.mask particle area using rectangle
% ָ����ת�Ƕȣ�˳ʱ��Ϊ������ʱ��Ϊ����
rotationAngle = 13; % ����Ը�����Ҫ������ת�Ƕ�

% �Բü����ͼƬ������ת
rotatedCroppedImageData = imrotate(croppedImageData, rotationAngle, 'bilinear', 'crop');

% ��ʾ��ת���ͼƬ
figure;  
imshow(rotatedCroppedImageData);
title('��ת���ͼƬ');

% ʹ��imrect�������þ��εĴ�С
disp('���ֶ����þ��εĴ�С...');
rect = imrect;
wait(rect);

% ��ȡ���ε�λ�úʹ�С
position = getPosition(rect);
x1 = round(position(1));
y1 = round(position(2));
width = round(position(3));
height = round(position(4));
x2 = x1 + width - 1;
y2 = y1 + height - 1;

% ����һ��������ģ
rectangleMask = zeros(size(rotatedCroppedImageData));
rectangleMask(y1:y2, x1:x2) = 1;

% ����ת��Ĳü�ͼƬӦ�þ�����ģ
maskedRegion = double(rotatedCroppedImageData) .* rectangleMask;

% ��ʾӦ���˾�����ģ��ͼ��
figure;
imshow(uint8(maskedRegion));
title('Ӧ���˾�����ģ��ͼ��');
%% ++.hanning window fft for mask area
% ��Ӧ����Բ����ģ��ͼ��Ӧ��Hanning��
hanningWindow = hann(size(maskedRegion, 1)) * hann(size(maskedRegion, 2))';
maskedRegion = maskedRegion .* hanningWindow;

% ����FFT�任
fftResult = fft2(maskedRegion);

% ��ʾFFT����ķ�����
figure;
imshow(log(abs(fftshift(fftResult)) + 1), []);
title('Բ����ģFFT�任�ķ�����');
%% 

