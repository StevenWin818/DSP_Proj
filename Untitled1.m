close all; clear; clc;

%% 线性调频(LFM)信号分析与滤波
% 信号参数
fc = 30e6;        % 载波频率 30MHz
B = 2e6;          % 带宽 2MHz
T = 300e-6;       % 脉宽 300us
K = B/T;          % 调频斜率

% 采样参数设置
fs = 10e6;        % 采样频率 10MHz (满足奈奎斯特采样定理)
t = -T/2:1/fs:T/2-1/fs; % 时间向量
N = length(t);    % 采样点数

%% 1. 生成LFM信号并进行时频分析
% 生成LFM信号 (实数形式)
s_real = cos(2*pi*(fc*t + 0.5*K*t.^2));

% 生成LFM信号 (复数形式，便于后续处理)
s_complex = exp(1j*2*pi*(fc*t + 0.5*K*t.^2));

% 绘制LFM信号的实部时域波形
figure('Name', 'LFM信号时域波形');
plot(t*1e6, real(s_complex));
xlabel('时间 (μs)');
ylabel('幅度');
title('LFM信号时域波形');
grid on;

% 计算并绘制频谱
freq = (-N/2:N/2-1)*(fs/N);
S = fftshift(fft(s_complex, N));

figure('Name', 'LFM信号频谱');
plot(freq/1e6, abs(S)/max(abs(S)));
xlabel('频率 (MHz)');
ylabel('归一化幅度');
title('LFM信号频谱');
xlim([-fs/2, fs/2]/1e6);
grid on;

% 时频分析 - 短时傅里叶变换 (STFT)
window_length = round(N/20);  % 窗长度
overlap = round(window_length*0.8);  % 重叠长度
nfft = 2^nextpow2(window_length*2);  % FFT点数

figure('Name', 'LFM信号时频图');
spectrogram(s_complex, hamming(window_length), overlap, nfft, fs, 'yaxis');
title('LFM信号时频图');

%% 2. 混频并添加噪声干扰
% 本地振荡器 (混频到中频2MHz)
f_if = 2e6;  % 中频
lo = exp(-1j*2*pi*(fc-f_if)*t);
s_if = s_complex .* lo;  % 混频后的信号

% 添加高斯白噪声
snr_dB = 10;  % 信噪比(dB)
s_if_noisy = awgn(s_if, snr_dB, 'measured');

% 绘制混频和加噪后的时域波形
figure('Name', '混频和加噪后信号');
subplot(2,1,1);
plot(t*1e6, real(s_if));
xlabel('时间 (μs)');
ylabel('幅度');
title('混频后信号的时域波形');
grid on;

subplot(2,1,2);
plot(t*1e6, real(s_if_noisy));
xlabel('时间 (μs)');
ylabel('幅度');
title(['混频后加噪信号的时域波形 (SNR = ', num2str(snr_dB), ' dB)']);
grid on;

% 频谱分析
S_if = fftshift(fft(s_if, N));
S_if_noisy = fftshift(fft(s_if_noisy, N));

figure('Name', '混频和加噪后频谱');
subplot(2,1,1);
plot(freq/1e6, abs(S_if)/max(abs(S_if)));
xlabel('频率 (MHz)');
ylabel('归一化幅度');
title('混频后信号的频谱');
grid on;
xlim([-5, 5]);  % 显示±5MHz范围

subplot(2,1,2);
plot(freq/1e6, abs(S_if_noisy)/max(abs(S_if_noisy)));
xlabel('频率 (MHz)');
ylabel('归一化幅度');
title(['混频后加噪信号的频谱 (SNR = ', num2str(snr_dB), ' dB)']);
grid on;
xlim([-fs/2e6, fs/2e6]);  % 显示可观测的频率范围(-5MHz到5MHz)

%% 3. 设计合适的滤波器
% 设计带通滤波器，覆盖LFM信号的频带
f_low = (f_if - B/2) / (fs/2);  % 归一化截止频率
f_high = (f_if + B/2) / (fs/2);  % 归一化截止频率

% 使用FIR滤波器（窗函数法）
filter_order = 128;  % 滤波器阶数
b_fir = fir1(filter_order, [f_low, f_high], 'bandpass', hamming(filter_order+1));
[h_fir, w_fir] = freqz(b_fir, 1, 1024);

% 绘制滤波器的幅度响应和相位响应
figure('Name', '带通滤波器响应');
subplot(2,1,1);
plot(w_fir/(2*pi)*fs/1e6, 20*log10(abs(h_fir)));
xlabel('频率 (MHz)');
ylabel('增益 (dB)');
title('带通滤波器幅度响应');
grid on;

subplot(2,1,2);
plot(w_fir/(2*pi)*fs/1e6, unwrap(angle(h_fir)));
xlabel('频率 (MHz)');
ylabel('相位 (弧度)');
title('带通滤波器相位响应');
grid on;

%% 4. 滤波处理并比较滤波前后效果
% 对带噪声的信号进行滤波
s_if_filtered = filter(b_fir, 1, s_if_noisy);

% 计算滤波前后的频谱
S_if_filtered = fftshift(fft(s_if_filtered, N));

% 在同一窗口中比较滤波前后的波形和频谱
figure('Name', '滤波效果比较');
% 时域波形比较
subplot(2,2,1);
plot(t*1e6, real(s_if));
hold on;
plot(t*1e6, real(s_if_noisy), 'r:');
xlabel('时间 (μs)');
ylabel('幅度');
title('原始信号和带噪声信号');
legend('原始信号', '带噪声信号');
grid on;

subplot(2,2,3);
plot(t*1e6, real(s_if));
hold on;
plot(t*1e6, real(s_if_filtered), 'g-.');
xlabel('时间 (μs)');
ylabel('幅度');
title('原始信号和滤波后信号');
legend('原始信号', '滤波后信号');
grid on;

% 频域比较
subplot(2,2,2);
plot(freq/1e6, abs(S_if)/max(abs(S_if)));
hold on;
plot(freq/1e6, abs(S_if_noisy)/max(abs(S_if_noisy)), 'r:');
xlabel('频率 (MHz)');
ylabel('归一化幅度');
title('原始信号和带噪声信号的频谱');
legend('原始信号', '带噪声信号');
grid on;
xlim([-5, 5]);  % 显示±5MHz范围

subplot(2,2,4);
plot(freq/1e6, abs(S_if)/max(abs(S_if)));
hold on;
plot(freq/1e6, abs(S_if_filtered)/max(abs(S_if_filtered)), 'g-.');
xlabel('频率 (MHz)');
ylabel('归一化幅度');
title('原始信号和滤波后信号的频谱');
legend('原始信号', '滤波后信号');
grid on;
xlim([-5, 5]);  % 显示±5MHz范围

% 计算滤波前后的信噪比改善
power_signal = sum(abs(s_if).^2)/N;
power_noise_before = sum(abs(s_if_noisy - s_if).^2)/N;
power_noise_after = sum(abs(s_if_filtered - s_if).^2)/N;

snr_before = 10*log10(power_signal/power_noise_before);
snr_after = 10*log10(power_signal/power_noise_after);
snr_improvement = snr_after - snr_before;

fprintf('滤波前信噪比: %.2f dB\n', snr_before);
fprintf('滤波后信噪比: %.2f dB\n', snr_after);
fprintf('信噪比改善: %.2f dB\n', snr_improvement);

% 匹配滤波分析（可选）
mf = conj(fliplr(s_complex));  % 匹配滤波器脉冲响应
y_mf_noisy = conv(s_if_noisy, mf);  % 匹配滤波输出
y_mf_filtered = conv(s_if_filtered, mf);  % 滤波后信号的匹配滤波输出

% 绘制匹配滤波输出
t_mf = (-N+1:N-1)/fs;  % 匹配滤波后的时间轴
figure('Name', '匹配滤波输出比较');
subplot(2,1,1);
plot(t_mf*1e6, abs(y_mf_noisy)/max(abs(y_mf_noisy)));
xlabel('时间 (μs)');
ylabel('归一化幅度');
title('带噪声信号的匹配滤波输出');
grid on;

subplot(2,1,2);
plot(t_mf*1e6, abs(y_mf_filtered)/max(abs(y_mf_filtered)));
xlabel('时间 (μs)');
ylabel('归一化幅度');
title('滤波后信号的匹配滤波输出');
grid on;