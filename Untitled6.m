%% 语音信号800Hz单频干扰去除 - 数字滤波器设计与分析
close all; clear; clc;

%% 1. 基本参数设置
fs = 8000;                % 采样频率 (Hz)
f_interfere = 800;        % 干扰频率 (Hz)
duration = 3;             % 信号持续时间 (秒)

% 生成时间向量
t = (0:1/fs:duration-1/fs)';

%% 2. 语音数据获取和准备
try
    % 尝试读取语音文件（如果存在）
    [speech, fs_orig] = audioread('speech_sample.wav');
    
    % 如果采样率与设定不同，需要重采样
    if fs_orig ~= fs
        speech = resample(speech, fs, fs_orig);
    end
    
    % 如果是立体声，转为单声道
    if size(speech, 2) > 1
        speech = mean(speech, 2);
    end
    
    % 裁剪或补长到所需长度
    if length(speech) > length(t)
        speech = speech(1:length(t));
    else
        speech = [speech; zeros(length(t)-length(speech), 1)];
    end
    
    % 归一化
    speech = speech / max(abs(speech));
    
catch
    % 如果没有语音文件，生成模拟语音信号（调制噪声）
    disp('找不到语音文件，使用模拟语音信号...');
    noise = randn(length(t), 1);
    
    % 低通滤波以模拟语音信号的频谱特性
    [b, a] = butter(6, 3000/(fs/2));
    speech = filter(b, a, noise);
    
    % 归一化
    speech = 0.8 * speech / max(abs(speech));
end

% 生成干扰信号
interference = 0.5 * sin(2*pi*f_interfere*t);

% 添加干扰到语音信号
noisy_speech = speech + interference;

%% 3. 滤波器设计
% 方法1: IIR陷波滤波器设计
Q_values = [10, 30, 50];  % 不同的品质因子
bw_values = f_interfere ./ Q_values;  % 带宽值

% 存储不同Q值的滤波器系数
b_notch = cell(length(Q_values), 1);
a_notch = cell(length(Q_values), 1);

for i = 1:length(Q_values)
    [b_notch{i}, a_notch{i}] = iirnotch(f_interfere/(fs/2), bw_values(i)/(fs/2));
end

% 方法2: FIR带阻滤波器设计（窗函数法）
filter_orders = [32, 64, 128];  % 不同的滤波器阶数
stop_width_values = [50, 100, 200];  % 不同的阻带宽度 (Hz)

% 存储不同参数的FIR滤波器系数
b_fir = cell(length(filter_orders), length(stop_width_values));

for i = 1:length(filter_orders)
    for j = 1:length(stop_width_values)
        % 设计带阻滤波器
        f_low = (f_interfere - stop_width_values(j)/2) / (fs/2);  % 归一化下截止频率
        f_high = (f_interfere + stop_width_values(j)/2) / (fs/2); % 归一化上截止频率
        
        % 确保频率在有效范围内
        f_low = max(0.001, min(0.999, f_low));
        f_high = max(0.001, min(0.999, f_high));
        
        b_fir{i,j} = fir1(filter_orders(i), [f_low f_high], 'stop', hamming(filter_orders(i)+1));
    end
end

%% 4. 滤波器频率响应分析
% 分析IIR陷波滤波器
figure('Name', 'IIR陷波滤波器频率响应');
[h1, w1] = freqz(b_notch{1}, a_notch{1}, 1024, fs);
[h2, w2] = freqz(b_notch{2}, a_notch{2}, 1024, fs);
[h3, w3] = freqz(b_notch{3}, a_notch{3}, 1024, fs);

subplot(2,1,1);
plot(w1, 20*log10(abs(h1)), 'b', w2, 20*log10(abs(h2)), 'r', w3, 20*log10(abs(h3)), 'g');
title('IIR陷波滤波器幅频响应');
xlabel('频率 (Hz)');
ylabel('增益 (dB)');
legend(['Q = ' num2str(Q_values(1))], ['Q = ' num2str(Q_values(2))], ['Q = ' num2str(Q_values(3))]);
grid on;
xlim([f_interfere-500 f_interfere+500]);

subplot(2,1,2);
plot(w1, unwrap(angle(h1))*180/pi, 'b', w2, unwrap(angle(h2))*180/pi, 'r', w3, unwrap(angle(h3))*180/pi, 'g');
title('IIR陷波滤波器相频响应');
xlabel('频率 (Hz)');
ylabel('相位 (度)');
grid on;
xlim([f_interfere-500 f_interfere+500]);

% 分析FIR带阻滤波器
figure('Name', 'FIR带阻滤波器频率响应');
[h_fir1, w_fir1] = freqz(b_fir{3,1}, 1, 1024, fs);
[h_fir2, w_fir2] = freqz(b_fir{3,2}, 1, 1024, fs);
[h_fir3, w_fir3] = freqz(b_fir{3,3}, 1, 1024, fs);

subplot(2,1,1);
plot(w_fir1, 20*log10(abs(h_fir1)), 'b', w_fir2, 20*log10(abs(h_fir2)), 'r', w_fir3, 20*log10(abs(h_fir3)), 'g');
title(['FIR带阻滤波器幅频响应 (阶数 = ' num2str(filter_orders(3)) ')']);
xlabel('频率 (Hz)');
ylabel('增益 (dB)');
legend(['阻带宽度 = ' num2str(stop_width_values(1)) ' Hz'], ...
       ['阻带宽度 = ' num2str(stop_width_values(2)) ' Hz'], ...
       ['阻带宽度 = ' num2str(stop_width_values(3)) ' Hz']);
grid on;
xlim([f_interfere-500 f_interfere+500]);

subplot(2,1,2);
plot(w_fir1, unwrap(angle(h_fir1))*180/pi, 'b', w_fir2, unwrap(angle(h_fir2))*180/pi, 'r', w_fir3, unwrap(angle(h_fir3))*180/pi, 'g');
title('FIR带阻滤波器相频响应');
xlabel('频率 (Hz)');
ylabel('相位 (度)');
grid on;
xlim([f_interfere-500 f_interfere+500]);

% 分析不同阶数FIR滤波器的影响
figure('Name', 'FIR带阻滤波器阶数影响');
[h_order1, w_order1] = freqz(b_fir{1,2}, 1, 1024, fs);
[h_order2, w_order2] = freqz(b_fir{2,2}, 1, 1024, fs);
[h_order3, w_order3] = freqz(b_fir{3,2}, 1, 1024, fs);

subplot(2,1,1);
plot(w_order1, 20*log10(abs(h_order1)), 'b', w_order2, 20*log10(abs(h_order2)), 'r', w_order3, 20*log10(abs(h_order3)), 'g');
title(['FIR带阻滤波器阶数对比 (阻带宽度 = ' num2str(stop_width_values(2)) ' Hz)']);
xlabel('频率 (Hz)');
ylabel('增益 (dB)');
legend(['阶数 = ' num2str(filter_orders(1))], ['阶数 = ' num2str(filter_orders(2))], ['阶数 = ' num2str(filter_orders(3))]);
grid on;
xlim([f_interfere-500 f_interfere+500]);

subplot(2,1,2);
plot(w_order1, unwrap(angle(h_order1))*180/pi, 'b', w_order2, unwrap(angle(h_order2))*180/pi, 'r', w_order3, unwrap(angle(h_order3))*180/pi, 'g');
title('FIR带阻滤波器阶数对相位的影响');
xlabel('频率 (Hz)');
ylabel('相位 (度)');
grid on;
xlim([f_interfere-500 f_interfere+500]);

%% 5. 滤波效果测试与对比
% 应用IIR陷波滤波器
filtered_iir = filter(b_notch{2}, a_notch{2}, noisy_speech);

% 应用FIR带阻滤波器
filtered_fir = filter(b_fir{3,2}, 1, noisy_speech);

% 计算时域信号的频谱
N = length(noisy_speech);
freq = (0:N/2) * fs / N;

Y_speech = fft(speech, N);
Y_speech_mag = abs(Y_speech(1:N/2+1));

Y_noisy = fft(noisy_speech, N);
Y_noisy_mag = abs(Y_noisy(1:N/2+1));

Y_filtered_iir = fft(filtered_iir, N);
Y_filtered_iir_mag = abs(Y_filtered_iir(1:N/2+1));

Y_filtered_fir = fft(filtered_fir, N);
Y_filtered_fir_mag = abs(Y_filtered_fir(1:N/2+1));

% 绘制频谱对比
figure('Name', '滤波前后的信号频谱对比');
subplot(2,1,1);
plot(freq, Y_speech_mag, 'b', freq, Y_noisy_mag, 'r');
title('原始语音与带干扰语音的频谱');
xlabel('频率 (Hz)');
ylabel('幅度');
legend('原始语音', '带干扰语音');
grid on;
xlim([0 fs/2]);
ylim([0 max(Y_noisy_mag)*1.1]);

% 标注干扰频率
hold on;
%h = line([f_interfere f_interfere], [0 max(Y_noisy_mag)*1.1], 'Color', 'k', 'LineStyle', '--');
text(f_interfere+10, max(Y_noisy_mag)*0.9, ['干扰频率: ' num2str(f_interfere) ' Hz'], 'FontSize', 8);
hold off;

subplot(2,1,2);
plot(freq, Y_speech_mag, 'b', freq, Y_filtered_iir_mag, 'g', freq, Y_filtered_fir_mag, 'r');
title('原始语音与滤波后语音的频谱');
xlabel('频率 (Hz)');
ylabel('幅度');
legend('原始语音', 'IIR陷波滤波', 'FIR带阻滤波');
grid on;
xlim([0 fs/2]);
ylim([0 max(Y_speech_mag)*1.1]);

% 标注干扰频率
hold on;
h = line([f_interfere f_interfere], [0 max(Y_speech_mag)*1.1], 'Color', 'k', 'LineStyle', '--');
text(f_interfere+10, max(Y_speech_mag)*0.9, ['干扰频率: ' num2str(f_interfere) ' Hz'], 'FontSize', 8);
hold off;

% 绘制时域信号对比
figure('Name', '时域信号对比');
subplot(4,1,1);
plot(t, speech);
title('原始语音');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;
xlim([0 min(3, duration)]);

subplot(4,1,2);
plot(t, noisy_speech);
title('带干扰语音');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;
xlim([0 min(3, duration)]);

subplot(4,1,3);
plot(t, filtered_iir);
title('IIR陷波滤波后语音');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;
xlim([0 min(3, duration)]);

subplot(4,1,4);
plot(t, filtered_fir);
title('FIR带阻滤波后语音');
xlabel('时间 (秒)');
ylabel('幅度');
grid on;
xlim([0 min(3, duration)]);

% 计算信噪比改善
power_signal = sum(speech.^2) / length(speech);
power_noise_before = sum((noisy_speech - speech).^2) / length(speech);
power_noise_iir = sum((filtered_iir - speech).^2) / length(speech);
power_noise_fir = sum((filtered_fir - speech).^2) / length(speech);

SNR_before = 10*log10(power_signal / power_noise_before);
SNR_iir = 10*log10(power_signal / power_noise_iir);
SNR_fir = 10*log10(power_signal / power_noise_fir);

fprintf('滤波前信噪比: %.2f dB\n', SNR_before);
fprintf('IIR陷波滤波后信噪比: %.2f dB (改善 %.2f dB)\n', SNR_iir, SNR_iir-SNR_before);
fprintf('FIR带阻滤波后信噪比: %.2f dB (改善 %.2f dB)\n', SNR_fir, SNR_fir-SNR_before);

%% 6. 播放音频比较效果（可选）
% 注释此部分如不想播放声音
disp('播放原始语音...');
sound(speech, fs);
pause(duration + 1);

disp('播放带干扰语音...');
sound(noisy_speech, fs);
pause(duration + 1);

disp('播放IIR陷波滤波后语音...');
sound(filtered_iir, fs);
pause(duration + 1);

disp('播放FIR带阻滤波后语音...');
sound(filtered_fir, fs);

% 保存处理后的音频（可选）
% audiowrite('original_speech.wav', speech, fs);
% audiowrite('noisy_speech.wav', noisy_speech, fs);
% audiowrite('filtered_iir.wav', filtered_iir, fs);
% audiowrite('filtered_fir.wav', filtered_fir, fs);