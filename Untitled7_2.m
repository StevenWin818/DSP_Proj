close all; clear; clc;

% --- 定义参数并调用函数 ---
audio_input_file = 'your_audio_file.wav'; % 替换为你的语音文件名，或保持此状态以使用自动生成的信号
noise_frequencies = [600, 1000, 1500];    % 要去除的多个干扰频率(Hz)

% 调用主函数
remove_multiple_sine_noise_dft(audio_input_file, noise_frequencies);

% -------------------------------
% 仅用DFT实现的去噪函数示例
function [clean_signal] = remove_multiple_sine_noise_dft(audio_file, noise_freqs)
    % 1. 读取或生成音频
    if exist(audio_file, 'file')
        [y, Fs] = audioread(audio_file);
        if size(y, 2) > 1
            y = y(:, 1); % 如果是立体声，只取第一个通道
        end
        disp(['已读取文件: ', audio_file]);
    else
        % 不存在音频文件，生成3秒1kHz正弦波模拟人声
        Fs = 16000; 
        t_gen = (0:Fs*3-1)'/Fs;
        y = 0.5 * sin(2*pi*1000*t_gen);
        disp('未找到语音文件，已自动生成3秒1kHz正弦波作为测试信号');
    end
    
    % 2. 合成带噪声信号
    t = (0:length(y)-1)'/Fs;
    noise = zeros(size(y));
    for i = 1:length(noise_freqs)
        noise = noise + 0.1 * sin(2*pi*noise_freqs(i)*t);
    end
    noisy_signal = y + noise;
    
    % 3. 对带噪声信号进行DFT
    N = length(noisy_signal);
    Y_dft = myDFT(noisy_signal);  % 自定义DFT
    
    % 4. 去除干扰频率
    df = Fs / N; % 频率分辨率
    Y_clean = Y_dft; 
    search_range = -3:3; % 在±3个频率点内搜索
    for k = 1:length(noise_freqs)
        % 估计干扰索引(正频率)
        freq_index_guess = round(noise_freqs(k)/df) + 1;
        if freq_index_guess < 1 || freq_index_guess > N
            continue;
        end
        
        % 在小范围内搜索峰值
        best_index = freq_index_guess;
        best_amp = abs(Y_dft(freq_index_guess));
        for s = search_range
            check_idx = freq_index_guess + s;
            if (check_idx >= 1 && check_idx <= N)
                amp_here = abs(Y_dft(check_idx));
                if amp_here > best_amp
                    best_amp = amp_here;
                    best_index = check_idx;
                end
            end
        end
        
        % 对找到峰值附近置零(正频率)
        bandwidth = 2; % 带宽
        lower_bound = max(best_index - bandwidth, 1);
        upper_bound = min(best_index + bandwidth, N);
        Y_clean(lower_bound:upper_bound) = 0;
        
        % 求对应的负频率索引
        neg_idx = N - best_index + 2;
        if neg_idx >= 1 && neg_idx <= N
            lower_bound_neg = max(neg_idx - bandwidth, 1);
            upper_bound_neg = min(neg_idx + bandwidth, N);
            Y_clean(lower_bound_neg:upper_bound_neg) = 0;
        end
    end
    
    % 5. 逆DFT获取去噪时域信号
    clean_signal_raw = real(myIDFT(Y_clean)); 

    % 6. 能量补偿，让去噪后信号能量与原信号保持接近
    energy_orig = sum(y.^2);
    energy_clean = sum(clean_signal_raw.^2);
    if energy_clean > 0
        scaling_factor = sqrt(energy_orig / energy_clean);
        clean_signal = clean_signal_raw * scaling_factor;
    else
        clean_signal = clean_signal_raw;
    end
    
    % 7. 幅度限制：防止去噪后信号截幅
    max_amp_orig = max(abs(y));
    max_amp_clean = max(abs(clean_signal));
    if max_amp_clean > max_amp_orig
        clean_signal = clean_signal * (max_amp_orig / max_amp_clean);
    end
    
    % 8. 绘图对比
    figure;
    subplot(3,1,1);
    plot(t, y);
    axis tight;
    title('原始语音信号');

    subplot(3,1,2);
    plot(t, noisy_signal);
    axis tight;
    title(['带多个正弦波干扰的信号 (', num2str(noise_freqs), 'Hz)']);

    subplot(3,1,3);
    plot(t, clean_signal);
    axis tight;
    title('去除多个正弦波干扰后的信号');
    
    % 频谱对比
    f = (0:N-1) * (Fs/N);
    figure;
    % 原始信号频谱
    Y_orig = myDFT(y);
    subplot(3,1,1);
    plot(f(1:floor(N/2)), abs(Y_orig(1:floor(N/2))));
    title('原始信号的频谱');
    
    % 带噪声频谱
    subplot(3,1,2);
    plot(f(1:floor(N/2)), abs(Y_dft(1:floor(N/2))));
    title('带噪声信号频谱');
    
    % 去噪后频谱
    subplot(3,1,3);
    plot(f(1:floor(N/2)), abs(Y_clean(1:floor(N/2))));
    title('去噪后信号频谱');
    
    % 播放对比
    disp('播放原始信号...');
    sound(y, Fs);
    pause(length(y)/Fs + 1);
    
    disp('播放带噪声信号...');
    sound(noisy_signal, Fs);
    pause(length(noisy_signal)/Fs + 1);
    
    disp('播放去噪后信号...');
    sound(clean_signal, Fs);
end

% -------------------------------
% 自定义DFT函数  
function X = myDFT(x)
    N = length(x);
    X = zeros(N, 1);
    for k = 1:N
        % 下面(k-1)和(n-1)是因为MATLAB索引从1开始
        for n = 1:N
            X(k) = X(k) + x(n)*exp(-1i*2*pi*(k-1)*(n-1)/N);
        end
    end
end

% -------------------------------
% 自定义IDFT函数
function x_rec = myIDFT(X)
    N = length(X);
    x_rec = zeros(N, 1);
    for n = 1:N
        for k = 1:N
            x_rec(n) = x_rec(n) + X(k)*exp(1i*2*pi*(k-1)*(n-1)/N);
        end
    end
    x_rec = x_rec / N; % IDFT要除以N
end