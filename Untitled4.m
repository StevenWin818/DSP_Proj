close all; clear; clc;

%% LFM信号调制与解调分析
% 参数定义
fc = 30e6;        % 载波频率 30MHz
B = 2e6;          % 带宽 2MHz
T = 300e-6;       % 脉宽 300us
K = B/T;          % 调频斜率

% 调制与解调参数
f_mod = 3e9;      % 调制载波频率 (3GHz)
f_demod = 2e9;    % 解调载波频率 (2GHz)

% 采样参数设置
% 注意：实际上为了完全表示3GHz载波需要更高采样率，这里我们使用相对低的采样率进行演示
% 实际应用中，可以使用复信号分析或降采样技术
fs = 10e9;        % 采样频率 10GHz (使高频载波可见)
dt = 1/fs;        % 采样间隔
t = -T/2:dt:T/2-dt; % 时间向量

%% 1. 产生LFM基带信号
% 生成LFM复信号
s_baseband = exp(1j*2*pi*(fc*t + 0.5*K*t.^2));

% 提取实部作为实信号分析
s_baseband_real = real(s_baseband);

% 频谱分析
N = length(t);
freq = (-N/2:N/2-1)*(fs/N); % 频率向量
S_baseband = fftshift(fft(s_baseband, N))/N;
S_baseband_real = fftshift(fft(s_baseband_real, N))/N;

% 绘制LFM基带信号的时域和频域图
figure('Name', 'LFM基带信号分析');
subplot(2,2,1);
plot(t*1e6, real(s_baseband));
title('LFM基带信号(复信号实部)');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

subplot(2,2,2);
plot(t*1e6, s_baseband_real);
title('LFM基带信号(实信号)');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

subplot(2,2,3);
plot(freq/1e6, abs(S_baseband));
title('LFM基带信号频谱(复信号)');
xlabel('频率 (MHz)');
ylabel('幅度');
grid on;
xlim([-100, 100]); % 限制频率轴显示范围

subplot(2,2,4);
plot(freq/1e6, abs(S_baseband_real));
title('LFM基带信号频谱(实信号)');
xlabel('频率 (MHz)');
ylabel('幅度');
grid on;
xlim([-100, 100]);

%% 2. 对基带信号进行3GHz载波调制
% 生成调制载波
carrier_mod = exp(1j*2*pi*f_mod*t);

% 调制LFM信号(使用复信号)
s_modulated = s_baseband .* carrier_mod;

% 转换为实信号
s_modulated_real = real(s_modulated);

% 频谱分析
S_modulated = fftshift(fft(s_modulated, N))/N;
S_modulated_real = fftshift(fft(s_modulated_real, N))/N;

% 绘制调制后信号的时域和频域特性
figure('Name', '3GHz调制信号分析');
subplot(2,2,1);
% 展示一小段调制信号以便观察高频载波
t_zoom_idx = 1:min(1000, length(t));
plot(t(t_zoom_idx)*1e9, real(s_modulated(t_zoom_idx)));
title('调制信号小段(复信号实部)');
xlabel('时间 (ns)');
ylabel('幅度');
grid on;

subplot(2,2,2);
plot(t(t_zoom_idx)*1e9, s_modulated_real(t_zoom_idx));
title('调制信号小段(实信号)');
xlabel('时间 (ns)');
ylabel('幅度');
grid on;

subplot(2,2,3);
plot(freq/1e9, abs(S_modulated));
title('调制信号频谱(复信号)');
xlabel('频率 (GHz)');
ylabel('幅度');
grid on;
xlim([-4, 4]); % 限制频率轴显示范围

subplot(2,2,4);
plot(freq/1e9, abs(S_modulated_real));
title('调制信号频谱(实信号)');
xlabel('频率 (GHz)');
ylabel('幅度');
grid on;
xlim([-4, 4]);

%% 3. 用2GHz载波解调信号
% 生成解调载波
carrier_demod = exp(-1j*2*pi*f_demod*t);

% 解调信号(使用复信号)
s_demodulated = s_modulated .* carrier_demod;

% 转换为实信号
s_demodulated_real = real(s_demodulated);

% 低通滤波以移除高频分量
% 设计低通滤波器
cutoff_freq = 100e6/(fs/2); % 截止频率100MHz
filter_order = 64;
b = fir1(filter_order, cutoff_freq, 'low');

% 滤波
s_demodulated_filtered = filter(b, 1, s_demodulated);
s_demodulated_real_filtered = filter(b, 1, s_demodulated_real);

% 频谱分析
S_demodulated = fftshift(fft(s_demodulated, N))/N;
S_demodulated_real = fftshift(fft(s_demodulated_real, N))/N;
S_demodulated_filtered = fftshift(fft(s_demodulated_filtered, N))/N;
S_demodulated_real_filtered = fftshift(fft(s_demodulated_real_filtered, N))/N;

% 绘制解调后信号的时域和频域特性
figure('Name', '2GHz解调信号分析');
subplot(2,2,1);
plot(t*1e6, real(s_demodulated_filtered));
title('滤波后解调信号(复信号实部)');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

subplot(2,2,2);
plot(t*1e6, s_demodulated_real_filtered);
title('滤波后解调信号(实信号)');
xlabel('时间 (μs)');
ylabel('幅度');
grid on;

subplot(2,2,3);
plot(freq/1e6, abs(S_demodulated));
title('解调信号频谱(滤波前)');
xlabel('频率 (MHz)');
ylabel('幅度');
grid on;
xlim([-2000, 2000]);

subplot(2,2,4);
plot(freq/1e6, abs(S_demodulated_filtered));
title('解调信号频谱(滤波后)');
xlabel('频率 (MHz)');
ylabel('幅度');
grid on;
xlim([-200, 200]);

%% 4. 分析调制与解调过程中的频率偏移
figure('Name', '调制解调过程频谱对比');
subplot(3,1,1);
plot(freq/1e6, abs(S_baseband));
title('基带LFM信号频谱');
xlabel('频率 (MHz)');
ylabel('幅度');
xlim([-100, 100]);
grid on;

subplot(3,1,2);
plot(freq/1e9, abs(S_modulated));
title('3GHz调制后信号频谱');
xlabel('频率 (GHz)');
ylabel('幅度');
xlim([-4, 4]);
grid on;

subplot(3,1,3);
plot(freq/1e6, abs(S_demodulated_filtered));
title('2GHz解调后滤波信号频谱');
xlabel('频率 (MHz)');
ylabel('幅度');
xlim([-1500, 1500]);
grid on;

% 显示频率偏移
freq_shift = f_mod - f_demod;
fprintf('调制频率与解调频率的差值: %.2f GHz\n', freq_shift/1e9);
fprintf('这导致解调信号相对于基带信号有 %.2f GHz 的频率偏移\n', freq_shift/1e9);