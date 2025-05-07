close all; clear; clc;

audio_file = 'your_audio_file.wav'; % 替换为你的语音文件名
noise_freq = 600; % 干扰频率

remove_sine_noise(audio_file, noise_freq);

function [clean_signal] = remove_sine_noise(audio_file, noise_freq)
% REMOVE_SINE_NOISE 使用DFT去除语音信号中的正弦波干扰
%   audio_file: 语音文件路径
%   noise_freq: 噪声频率(Hz)
%   clean_signal: 去除噪声后的信号

% 检查文件是否存在
if exist(audio_file, 'file')
    [y, Fs] = audioread(audio_file);
    if size(y, 2) > 1
        y = y(:, 1); % 如果是立体声，只取第一个通道
    end
    disp(['已读取文件: ', audio_file]);
else
    % 文件不存在，生成3秒1kHz正弦波模拟语音
    Fs = 16000; % 采样率
    t = (0:Fs*3-1)'/Fs;
    y = 0.5 * sin(2*pi*1000*t);
    disp('未找到语音文件，已自动生成3秒1kHz正弦波作为测试信号。');
end
t = (0:length(y)-1)'/Fs; % 时间向量

% 生成正弦波噪声
noise = 0.1 * sin(2*pi*noise_freq*t); % 幅度0.1的正弦波噪声

% 添加噪声到原始信号
noisy_signal = y + noise;

% 对带噪声信号进行FFT
N = length(noisy_signal);
Y = fft(noisy_signal);

% 计算频率分辨率和频率向量
df = Fs/N;
f = (0:N-1)*df;

% 找到噪声频率对应的索引
freq_index = round(noise_freq/df) + 1;
bandwidth = 5; % 设置带宽，去除噪声频率附近的成分

% 在频域中去除噪声
Y_clean = Y;
% 去除正频率部分
lower_bound = max(1, freq_index - bandwidth);
upper_bound = min(N, freq_index + bandwidth);
Y_clean(lower_bound:upper_bound) = 0;

% 去除负频率部分
neg_freq_index = N - freq_index + 2;
if neg_freq_index <= N
    lower_bound = max(1, neg_freq_index - bandwidth);
    upper_bound = min(N, neg_freq_index + bandwidth);
    Y_clean(lower_bound:upper_bound) = 0;
end

% 反FFT获得去噪后的信号
clean_signal = real(ifft(Y_clean));

% 绘制原始信号、带噪声信号和去噪后的信号的波形
figure;
subplot(3,1,1);
plot(t, y);
axis tight; % 自动适应幅度
title('原始语音信号');
xlabel('时间 (秒)');
ylabel('幅度');

subplot(3,1,2);
plot(t, noisy_signal);
axis tight;
title(['带', num2str(noise_freq), 'Hz正弦波干扰的信号']);
xlabel('时间 (秒)');
ylabel('幅度');

subplot(3,1,3);
plot(t, clean_signal);
axis tight;
title(['去除', num2str(noise_freq), 'Hz干扰后的信号']);
xlabel('时间 (秒)');
ylabel('幅度');

% 绘制频谱图
figure;
subplot(3,1,1);
Yy = fft(y); % 先计算fft
plot(f(1:N/2), abs(Yy(1:N/2)));
title('原始信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅度');

subplot(3,1,2);
plot(f(1:N/2), abs(Y(1:N/2)));
title(['带', num2str(noise_freq), 'Hz干扰的信号频谱']);
xlabel('频率 (Hz)');
ylabel('幅度');

subplot(3,1,3);
plot(f(1:N/2), abs(Y_clean(1:N/2)));
title(['去除', num2str(noise_freq), 'Hz干扰后的信号频谱']);
xlabel('频率 (Hz)');
ylabel('幅度');

% 播放声音进行比较
disp('播放原始语音...');
sound(y, Fs);
pause(length(y)/Fs + 1);

disp('播放带噪声语音...');
sound(noisy_signal, Fs);
pause(length(noisy_signal)/Fs + 1);

disp('播放去噪后语音...');
sound(clean_signal, Fs);
end