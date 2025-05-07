% 信号处理实验：DFT分析和滤波器设计
% 通过FFT实现谱分解和设计滤波器分离信号
close all; clear; clc;

%% 1. 生成复合正弦信号并计算DFT
fs = 2000;                % 采样率 2000Hz
t_duration = 1;           % 信号持续时间1秒
t = 0:1/fs:t_duration-1/fs; % 时间向量
N = length(t);            % 采样点数

% 生成三个频率的正弦信号
f1 = 120;   % 第一个正弦信号频率 120Hz
f2 = 350;   % 第二个正弦信号频率 350Hz
f3 = 700;   % 第三个正弦信号频率 700Hz

% 分别生成三个正弦信号并合成复合信号
x1 = sin(2*pi*f1*t);
x2 = sin(2*pi*f2*t);
x3 = sin(2*pi*f3*t);
x = x1 + x2 + x3;

% 绘制复合信号的时域波形
figure;
subplot(4,1,1);
plot(t, x1);
title('120Hz正弦信号');
xlabel('时间 (秒)'); ylabel('幅度');

subplot(4,1,2);
plot(t, x2);
title('350Hz正弦信号');
xlabel('时间 (秒)'); ylabel('幅度');

subplot(4,1,3);
plot(t, x3);
title('700Hz正弦信号');
xlabel('时间 (秒)'); ylabel('幅度');

subplot(4,1,4);
plot(t, x);
title('复合信号');
xlabel('时间 (秒)'); ylabel('幅度');

% 计算DFT并绘制幅频特性
X = fft(x);         % 计算复合信号的DFT
X_mag = abs(X)/N;   % 幅度谱，除以N进行归一化
X_mag = X_mag(1:N/2+1); % 取一半频谱（由于共轭对称性）
X_mag(2:end-1) = 2*X_mag(2:end-1); % 调整幅度

% 计算频率轴
f = fs*(0:(N/2))/N;

% 绘制幅频特性
figure;
plot(f, X_mag);
title('复合信号的幅频特性');
xlabel('频率 (Hz)'); ylabel('幅度');
grid on;
axis([0 1000 0 1.5]); % 限制x轴范围便于观察

%% 2. 添加随机噪声并分析时频特性
% 生成随机噪声
noise_level = 0.5;
noise = noise_level * randn(size(t));
x_noisy = x + noise;

% 计算带噪声信号的DFT
X_noisy = fft(x_noisy);
X_noisy_mag = abs(X_noisy)/N;
X_noisy_mag = X_noisy_mag(1:N/2+1);
X_noisy_mag(2:end-1) = 2*X_noisy_mag(2:end-1);

% 绘制带噪声信号的时域和频域图
figure;
subplot(2,2,1);
plot(t, x);
title('原始复合信号');
xlabel('时间 (秒)'); ylabel('幅度');

subplot(2,2,2);
plot(t, x_noisy);
title('带噪声的复合信号');
xlabel('时间 (秒)'); ylabel('幅度');

subplot(2,2,3);
plot(f, X_mag);
title('原始信号的幅频特性');
xlabel('频率 (Hz)'); ylabel('幅度');
axis([0 1000 0 1.5]);

subplot(2,2,4);
plot(f, X_noisy_mag);
title('带噪声信号的幅频特性');
xlabel('频率 (Hz)'); ylabel('幅度');
axis([0 1000 0 1.5]);

% 时频分析
figure;
spectrogram(x_noisy, hamming(256), 128, 1024, fs, 'yaxis');
title('带噪声信号的时频图');

%% 3. 设计滤波器并分离信号
% 方法一：使用FIR滤波器（窗函数法）
order = 100;  % 滤波器阶数

% 120Hz带通滤波器
f1_low = 100/fs;   % 通带下限
f1_high = 140/fs;  % 通带上限
b1_fir = fir1(order, [f1_low f1_high], 'bandpass');
[h1_fir, w1_fir] = freqz(b1_fir, 1, 1024);

% 350Hz带通滤波器
f2_low = 330/fs;
f2_high = 370/fs;
b2_fir = fir1(order, [f2_low f2_high], 'bandpass');
[h2_fir, w2_fir] = freqz(b2_fir, 1, 1024);

% 700Hz带通滤波器
f3_low = 680/fs;
f3_high = 720/fs;
b3_fir = fir1(order, [f3_low f3_high], 'bandpass');
[h3_fir, w3_fir] = freqz(b3_fir, 1, 1024);

% 方法二：使用IIR滤波器（巴特沃斯滤波器）
order_iir = 4;
[b1_iir, a1_iir] = butter(order_iir, [f1_low f1_high], 'bandpass');
[h1_iir, w1_iir] = freqz(b1_iir, a1_iir, 1024);

[b2_iir, a2_iir] = butter(order_iir, [f2_low f2_high], 'bandpass');
[h2_iir, w2_iir] = freqz(b2_iir, a2_iir, 1024);

[b3_iir, a3_iir] = butter(order_iir, [f3_low f3_high], 'bandpass');
[h3_iir, w3_iir] = freqz(b3_iir, a3_iir, 1024);

%% 4. 显示滤波器的幅度响应和相位响应
% FIR滤波器响应
figure;
subplot(3,2,1);
plot(w1_fir*fs/(2*pi), abs(h1_fir));
title('120Hz FIR带通滤波器的幅频响应');
xlabel('频率 (Hz)'); ylabel('增益');

subplot(3,2,2);
plot(w1_fir*fs/(2*pi), unwrap(angle(h1_fir)));
title('120Hz FIR带通滤波器的相位响应');
xlabel('频率 (Hz)'); ylabel('相位 (弧度)');

subplot(3,2,3);
plot(w2_fir*fs/(2*pi), abs(h2_fir));
title('350Hz FIR带通滤波器的幅频响应');
xlabel('频率 (Hz)'); ylabel('增益');

subplot(3,2,4);
plot(w2_fir*fs/(2*pi), unwrap(angle(h2_fir)));
title('350Hz FIR带通滤波器的相位响应');
xlabel('频率 (Hz)'); ylabel('相位 (弧度)');

subplot(3,2,5);
plot(w3_fir*fs/(2*pi), abs(h3_fir));
title('700Hz FIR带通滤波器的幅频响应');
xlabel('频率 (Hz)'); ylabel('增益');

subplot(3,2,6);
plot(w3_fir*fs/(2*pi), unwrap(angle(h3_fir)));
title('700Hz FIR带通滤波器的相位响应');
xlabel('频率 (Hz)'); ylabel('相位 (弧度)');

% IIR滤波器响应
figure;
subplot(3,2,1);
plot(w1_iir*fs/(2*pi), abs(h1_iir));
title('120Hz IIR带通滤波器的幅频响应');
xlabel('频率 (Hz)'); ylabel('增益');

subplot(3,2,2);
plot(w1_iir*fs/(2*pi), unwrap(angle(h1_iir)));
title('120Hz IIR带通滤波器的相位响应');
xlabel('频率 (Hz)'); ylabel('相位 (弧度)');

subplot(3,2,3);
plot(w2_iir*fs/(2*pi), abs(h2_iir));
title('350Hz IIR带通滤波器的幅频响应');
xlabel('频率 (Hz)'); ylabel('增益');

subplot(3,2,4);
plot(w2_iir*fs/(2*pi), unwrap(angle(h2_iir)));
title('350Hz IIR带通滤波器的相位响应');
xlabel('频率 (Hz)'); ylabel('相位 (弧度)');

subplot(3,2,5);
plot(w3_iir*fs/(2*pi), abs(h3_iir));
title('700Hz IIR带通滤波器的幅频响应');
xlabel('频率 (Hz)'); ylabel('增益');

subplot(3,2,6);
plot(w3_iir*fs/(2*pi), unwrap(angle(h3_iir)));
title('700Hz IIR带通滤波器的相位响应');
xlabel('频率 (Hz)'); ylabel('相位 (弧度)');

%% 5. 滤波并比较效果
% 使用FIR滤波器滤波
y1_fir = filter(b1_fir, 1, x_noisy);
y2_fir = filter(b2_fir, 1, x_noisy);
y3_fir = filter(b3_fir, 1, x_noisy);

% 使用IIR滤波器滤波
y1_iir = filter(b1_iir, a1_iir, x_noisy);
y2_iir = filter(b2_iir, a2_iir, x_noisy);
y3_iir = filter(b3_iir, a3_iir, x_noisy);

% 绘制FIR滤波结果
figure;
subplot(3,1,1);
plot(t, y1_fir, t, x1, 'r--');
title('FIR滤波器提取的120Hz信号');
xlabel('时间 (秒)'); ylabel('幅度');
legend('滤波后', '原始信号');

subplot(3,1,2);
plot(t, y2_fir, t, x2, 'r--');
title('FIR滤波器提取的350Hz信号');
xlabel('时间 (秒)'); ylabel('幅度');
legend('滤波后', '原始信号');

subplot(3,1,3);
plot(t, y3_fir, t, x3, 'r--');
title('FIR滤波器提取的700Hz信号');
xlabel('时间 (秒)'); ylabel('幅度');
legend('滤波后', '原始信号');

% 绘制IIR滤波结果
figure;
subplot(3,1,1);
plot(t, y1_iir, t, x1, 'r--');
title('IIR滤波器提取的120Hz信号');
xlabel('时间 (秒)'); ylabel('幅度');
legend('滤波后', '原始信号');

subplot(3,1,2);
plot(t, y2_iir, t, x2, 'r--');
title('IIR滤波器提取的350Hz信号');
xlabel('时间 (秒)'); ylabel('幅度');
legend('滤波后', '原始信号');

subplot(3,1,3);
plot(t, y3_iir, t, x3, 'r--');
title('IIR滤波器提取的700Hz信号');
xlabel('时间 (秒)'); ylabel('幅度');
legend('滤波后', '原始信号');

% 计算滤波器性能指标
mse_fir = [mean((y1_fir - x1).^2), mean((y2_fir - x2).^2), mean((y3_fir - x3).^2)];
mse_iir = [mean((y1_iir - x1).^2), mean((y2_iir - x2).^2), mean((y3_iir - x3).^2)];

fprintf('FIR滤波器的均方误差：\n');
fprintf('120Hz: %.6f\n', mse_fir(1));
fprintf('350Hz: %.6f\n', mse_fir(2));
fprintf('700Hz: %.6f\n', mse_fir(3));

fprintf('\nIIR滤波器的均方误差：\n');
fprintf('120Hz: %.6f\n', mse_iir(1));
fprintf('350Hz: %.6f\n', mse_iir(2));
fprintf('700Hz: %.6f\n', mse_iir(3));