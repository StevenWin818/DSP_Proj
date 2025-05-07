%% 音频信号的采样及处理研究
close all; clear; clc;

%% 1. 以44kHz的采样频率采集音频信号
% 读取音频文件或生成测试信号
try
    % 尝试读取音频文件
    [y_orig, Fs_orig] = audioread('audio_sample.wav');
    
    % 如果是立体声，转换为单声道
    if size(y_orig, 2) > 1
        y_orig = mean(y_orig, 2);
    end
    
    % 如果原始采样率不是44kHz，重采样
    if Fs_orig ~= 44000
        y_44k = resample(y_orig, 44000, Fs_orig);
        Fs_44k = 44000;
    else
        y_44k = y_orig;
        Fs_44k = Fs_orig;
    end
    
    disp('成功读取音频文件');
catch
    % 如果没有找到音频文件，生成一个测试信号
    Fs_44k = 44000;
    duration = 3; % 3秒的信号
    t = 0:1/Fs_44k:(duration-1/Fs_44k);
    
    % 生成包含多个频率成分的测试信号
    y_44k = 0.5*sin(2*pi*440*t)' + 0.3*sin(2*pi*1000*t)' + 0.2*sin(2*pi*5000*t)' + 0.1*sin(2*pi*10000*t)';
    
    disp('使用生成的测试信号，包含440Hz, 1kHz, 5kHz和10kHz的频率成分');
end

% 绘制原始信号的波形和频谱
figure('Name', '原始音频信号（44kHz采样）');
subplot(2,1,1);
t_44k = (0:length(y_44k)-1)/Fs_44k;
plot(t_44k, y_44k);
title('原始音频信号波形（44kHz采样）');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;

subplot(2,1,2);
[pxx_44k, f_44k] = pwelch(y_44k, hamming(1024), 512, 1024, Fs_44k);
plot(f_44k, 10*log10(pxx_44k));
title('原始音频信号频谱（44kHz采样）');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
grid on;
xlim([0 Fs_44k/2]);

% 播放原始音频
disp('播放原始44kHz采样的音频...');
sound(y_44k, Fs_44k);
pause(length(y_44k)/Fs_44k + 0.5);

%% 2. 以4kHz的采样频率采集同一段音频，体会混叠现象
% 重采样到4kHz
Fs_4k = 4000;
y_4k = resample(y_44k, Fs_4k, Fs_44k);

% 绘制4kHz采样信号的波形和频谱
figure('Name', '重采样音频信号（4kHz采样）');
subplot(2,1,1);
t_4k = (0:length(y_4k)-1)/Fs_4k;
plot(t_4k, y_4k);
title('重采样音频信号波形（4kHz采样）');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;

subplot(2,1,2);
[pxx_4k, f_4k] = pwelch(y_4k, hamming(1024), 512, 1024, Fs_4k);
plot(f_4k, 10*log10(pxx_4k));
title('重采样音频信号频谱（4kHz采样）');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
grid on;
xlim([0 Fs_4k/2]);

% 播放4kHz采样的音频
disp('播放4kHz采样的音频（可能有混叠失真）...');
sound(y_4k, Fs_4k);
pause(length(y_4k)/Fs_4k + 0.5);

% 为了更直观地比较混叠现象，将4kHz采样的信号重采样回44kHz
y_4k_to_44k = resample(y_4k, Fs_44k, Fs_4k);
% 如果长度不一致，调整到相同长度
min_len = min(length(y_44k), length(y_4k_to_44k));
y_44k_cut = y_44k(1:min_len);
y_4k_to_44k_cut = y_4k_to_44k(1:min_len);

% 比较原始44kHz和从4kHz重采样回来的信号频谱
figure('Name', '混叠现象演示');
[pxx_44k_cut, f_44k_cut] = pwelch(y_44k_cut, hamming(1024), 512, 1024, Fs_44k);
[pxx_4k_to_44k, f_4k_to_44k] = pwelch(y_4k_to_44k_cut, hamming(1024), 512, 1024, Fs_44k);

plot(f_44k_cut, 10*log10(pxx_44k_cut), 'b', f_4k_to_44k, 10*log10(pxx_4k_to_44k), 'r');
title('混叠现象对比');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
legend('原始44kHz信号', '4kHz采样后重建的信号');
grid on;
xlim([0 20000]);

% 播放从4kHz重采样回44kHz的音频
disp('播放4kHz采样后重建回44kHz的音频...');
sound(y_4k_to_44k, Fs_44k);
pause(length(y_4k_to_44k)/Fs_44k + 0.5);

%% 3. 添加噪声并设计滤波器消除
% 添加不同类型的噪声
% 1. 白噪声
noise_level = 0.05;
y_noisy_white = y_44k + noise_level * randn(size(y_44k));

% 2. 高频噪声（模拟电器噪声）
t_noise = (0:length(y_44k)-1)/Fs_44k;
high_freq_noise = 0.1 * sin(2*pi*15000*t_noise)';
y_noisy_high = y_44k + high_freq_noise;

% 绘制带噪声信号的频谱
figure('Name', '带噪声信号频谱');

subplot(2,1,1);
[pxx_white, f_white] = pwelch(y_noisy_white, hamming(1024), 512, 1024, Fs_44k);
plot(f_white, 10*log10(pxx_white));
title('带白噪声的信号频谱');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
grid on;
xlim([0 Fs_44k/2]);

subplot(2,1,2);
[pxx_high, f_high] = pwelch(y_noisy_high, hamming(1024), 512, 1024, Fs_44k);
plot(f_high, 10*log10(pxx_high));
title('带高频噪声的信号频谱');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
grid on;
xlim([0 Fs_44k/2]);

% 设计滤波器消除噪声
% 1. 设计低通滤波器消除高频噪声
cutoff_freq = 10000/(Fs_44k/2);  % 10kHz截止频率
filter_order = 50;
b_lowpass = fir1(filter_order, cutoff_freq, 'low', hamming(filter_order+1));
y_filtered_high = filter(b_lowpass, 1, y_noisy_high);

% 2. 对于白噪声，使用带通滤波器保留人耳听觉范围内的频率
low_cutoff = 100/(Fs_44k/2);  % 100Hz低端截止
high_cutoff = 8000/(Fs_44k/2); % 8kHz高端截止
b_bandpass = fir1(filter_order, [low_cutoff high_cutoff], 'bandpass', hamming(filter_order+1));
y_filtered_white = filter(b_bandpass, 1, y_noisy_white);

% 绘制滤波器频率响应
figure('Name', '滤波器频率响应');

subplot(2,1,1);
[h_lowpass, w_lowpass] = freqz(b_lowpass, 1, 1024, Fs_44k);
plot(w_lowpass, 20*log10(abs(h_lowpass)));
title('低通滤波器频率响应');
xlabel('频率 (Hz)');
ylabel('幅度 (dB)');
grid on;

subplot(2,1,2);
[h_bandpass, w_bandpass] = freqz(b_bandpass, 1, 1024, Fs_44k);
plot(w_bandpass, 20*log10(abs(h_bandpass)));
title('带通滤波器频率响应');
xlabel('频率 (Hz)');
ylabel('幅度 (dB)');
grid on;

% 绘制滤波前后的频谱对比
figure('Name', '滤波前后频谱对比');

% 高频噪声滤波对比
subplot(2,2,1);
plot(f_high, 10*log10(pxx_high));
title('带高频噪声的信号频谱');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
grid on;
xlim([0 Fs_44k/2]);

subplot(2,2,2);
[pxx_filtered_high, f_filtered_high] = pwelch(y_filtered_high, hamming(1024), 512, 1024, Fs_44k);
plot(f_filtered_high, 10*log10(pxx_filtered_high));
title('低通滤波后的信号频谱');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
grid on;
xlim([0 Fs_44k/2]);

% 白噪声滤波对比
subplot(2,2,3);
plot(f_white, 10*log10(pxx_white));
title('带白噪声的信号频谱');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
grid on;
xlim([0 Fs_44k/2]);

subplot(2,2,4);
[pxx_filtered_white, f_filtered_white] = pwelch(y_filtered_white, hamming(1024), 512, 1024, Fs_44k);
plot(f_filtered_white, 10*log10(pxx_filtered_white));
title('带通滤波后的信号频谱');
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
grid on;
xlim([0 Fs_44k/2]);

% 时域波形对比
figure('Name', '滤波前后时域波形对比');
t_short = t_44k(1:min(10000, length(t_44k)));  % 只显示一小段以便观察细节

subplot(2,2,1);
plot(t_short, y_noisy_high(1:length(t_short)));
title('带高频噪声的时域波形');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;

subplot(2,2,2);
plot(t_short, y_filtered_high(1:length(t_short)));
title('低通滤波后的时域波形');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;

subplot(2,2,3);
plot(t_short, y_noisy_white(1:length(t_short)));
title('带白噪声的时域波形');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;

subplot(2,2,4);
plot(t_short, y_filtered_white(1:length(t_short)));
title('带通滤波后的时域波形');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;

% 计算信噪比改善
SNR_before_high = snr(y_44k, y_noisy_high - y_44k);
SNR_after_high = snr(y_44k, y_filtered_high - y_44k);
SNR_before_white = snr(y_44k, y_noisy_white - y_44k);
SNR_after_white = snr(y_44k, y_filtered_white - y_44k);

fprintf('高频噪声滤波前SNR: %.2f dB\n', SNR_before_high);
fprintf('高频噪声滤波后SNR: %.2f dB\n', SNR_after_high);
fprintf('白噪声滤波前SNR: %.2f dB\n', SNR_before_white);
fprintf('白噪声滤波后SNR: %.2f dB\n', SNR_after_white);

% 播放滤波前后的音频
disp('播放带高频噪声的音频...');
sound(y_noisy_high, Fs_44k);
pause(length(y_noisy_high)/Fs_44k + 0.5);

disp('播放低通滤波后的音频...');
sound(y_filtered_high, Fs_44k);
pause(length(y_filtered_high)/Fs_44k + 0.5);

disp('播放带白噪声的音频...');
sound(y_noisy_white, Fs_44k);
pause(length(y_noisy_white)/Fs_44k + 0.5);

disp('播放带通滤波后的音频...');
sound(y_filtered_white, Fs_44k);