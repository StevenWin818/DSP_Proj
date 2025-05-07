%% 线性卷积与圆周卷积关系、重叠相加法研究
close all; clear; clc;

%% 1. 产生两个离散时间信号并验证线性卷积和圆周卷积的关系
% 生成300点的信号x[n]和40点的信号h[n]
N1 = 300;
N2 = 40;

% 生成信号（使用随机数，也可以使用其他有意义的信号）
x = randn(N1, 1);       % 长信号, 列向量
h = fir1(N2-1, 0.25)';  % 短信号（使用低通滤波器作为冲激响应）, 转置为列向量以修正vertcat错误

% 计算线性卷积（使用MATLAB内置函数）
y_linear = conv(x, h);

% 计算圆周卷积
% 注意：直接用FFT乘积计算圆周卷积时，两个信号长度必须相同
% 情况1：直接计算圆周卷积（长度不足）
N_circ1 = N1;  % 使用长信号长度
x_padded1 = x;  % x已经是N1长度
h_padded1 = [h; zeros(N_circ1-N2, 1)];  % 对h进行零填充 (h是列向量，所以vertcat正确)
X1 = fft(x_padded1);
H1 = fft(h_padded1);
Y1 = X1 .* H1;
y_circ1 = real(ifft(Y1));  % 取实部以处理数值误差

% 情况2：计算足够长度的圆周卷积（等于线性卷积）
N_circ2 = N1 + N2 - 1;  % 线性卷积结果长度
x_padded2 = [x; zeros(N_circ2-N1, 1)];  % 对x进行零填充到足够长度
h_padded2 = [h; zeros(N_circ2-N2, 1)];  % 对h进行零填充到足够长度 (h是列向量)
X2 = fft(x_padded2);
H2 = fft(h_padded2);
Y2 = X2 .* H2;
y_circ2 = real(ifft(Y2));  % 取实部以处理数值误差

% 绘图比较线性卷积与圆周卷积
figure;
subplot(3,1,1);
stem(0:length(y_linear)-1, real(y_linear), 'filled');
title('线性卷积 y = x * h');
xlabel('样本点'); ylabel('幅度');
grid on;

subplot(3,1,2);
stem(0:length(y_circ1)-1, real(y_circ1), 'filled');
title(['长度为', num2str(N_circ1), '的圆周卷积（不足长度）']);
xlabel('样本点'); ylabel('幅度');
grid on;

subplot(3,1,3);
stem(0:length(y_circ2)-1, real(y_circ2), 'filled');
title(['长度为', num2str(N_circ2), '的圆周卷积（足够长度）']);
xlabel('样本点'); ylabel('幅度');
grid on;

% 计算线性卷积与圆周卷积的误差
error1 = max(abs(y_linear(1:N1) - y_circ1(1:N1))); % 比较相同长度部分
error2 = max(abs(y_linear - y_circ2));

fprintf('线性卷积与不足长度圆周卷积的最大误差 (比较前N1点): %.10f\n', error1);
fprintf('线性卷积与足够长度圆周卷积的最大误差: %.10f\n', error2);

%% 2. 编程实现重叠相加法
% 自定义函数实现重叠相加法
tic;
y_overlap = overlapAdd(x, h, 64); % 分段长度设为64
time_overlap = toc;

% 验证重叠相加法的正确性
error_overlap = max(abs(y_linear - y_overlap));
fprintf('重叠相加法与直接线性卷积的最大误差: %.10f\n', error_overlap);

%% 3. 比较计算量
% 直接计算线性卷积的时间
tic;
y_direct = conv(x, h);
time_direct = toc;

% 不同分段长度下重叠相加法的性能
segment_lengths = [32, 64, 128, 256, 512]; % 确保 L > N2 (40) for overlapAdd
valid_segment_lengths = segment_lengths(segment_lengths > N2);
if length(valid_segment_lengths) < length(segment_lengths)
    fprintf('警告: 部分分段长度 <= 滤波器长度，已从比较中移除。\n');
    segment_lengths = valid_segment_lengths;
end

times = zeros(size(segment_lengths));
for i = 1:length(segment_lengths)
    tic;
    overlapAdd(x, h, segment_lengths(i));
    times(i) = toc;
end

% 计算理论复杂度
N = N1;
M = N2;
direct_ops = N * M;  % 直接线性卷积的乘加次数
fft_ops = zeros(size(segment_lengths));

for i = 1:length(segment_lengths)
    L_seg = segment_lengths(i);
    P_fft = L_seg + M - 1;  % FFT长度
    % num_segments_calc = ceil(N / (L_seg - M + 1)); % 与overlapAdd函数中一致
    % 修正：重叠相加法，每段输入L，输出L+M-1，输出叠加时每隔L个点开始叠加
    num_segments_calc = ceil(N / L_seg); % 每个输入段长L
    % 每段FFT的复杂度(FFT+IFFT+乘法)
    % 实际操作中，通常用复数乘法估算，一个复数乘法约等于 3-4个实数乘法+若干加法
    % FFT/IFFT O(P_fft log P_fft)
    fft_ops(i) = num_segments_calc * (2 * 2.5 * P_fft * log2(P_fft) + P_fft); % 假设 FFT/IFFT 为 2.5*P*log2(P) 实数操作, P_fft复数乘法
end

% 绘制比较图表
figure;
subplot(2,1,1);
if ~isempty(times) % Ensure times is not empty before plotting
    bar_data = [time_direct, times];
    bar_labels = ['Direct', arrayfun(@(x) ['L=', num2str(x)], segment_lengths, 'UniformOutput', false)];
    bar(bar_data);
    set(gca, 'XTick', 1:length(bar_data));
    set(gca, 'XTickLabel', bar_labels);
else
    bar(time_direct); % Only direct if no valid segment lengths
    set(gca, 'XTick', 1);
    set(gca, 'XTickLabel', {'Direct'});
    disp('没有有效的分段长度用于重叠相加法性能比较。');
end
title('计算时间比较');
ylabel('时间 (秒)');
grid on;

subplot(2,1,2);
if ~isempty(fft_ops)
    ops_data = [direct_ops, fft_ops] / 1e6;
    ops_labels = ['Direct', arrayfun(@(x) ['L=', num2str(x)], segment_lengths, 'UniformOutput', false)];
    bar(ops_data);
    set(gca, 'XTick', 1:length(ops_data));
    set(gca, 'XTickLabel', ops_labels);
else
    bar(direct_ops / 1e6);
    set(gca, 'XTick', 1);
    set(gca, 'XTickLabel', {'Direct'});
end
title('理论计算量比较 (近似)');
ylabel('操作次数 (百万)');
grid on;

fprintf('直接计算线性卷积时间: %.6f 秒\n', time_direct);
for i = 1:length(segment_lengths)
    fprintf('重叠相加法(分段长度=%d)时间: %.6f 秒\n', segment_lengths(i), times(i));
end

%% 自定义重叠相加法函数
function y = overlapAdd(x, h, L)
    % 重叠相加法计算线性卷积
    % x: 输入信号 (列向量)
    % h: 滤波器冲激响应 (列向量)
    % L: 分段长度 (x的每个分段的长度)
    
    N = length(x);    % 输入信号长度
    M = length(h);    % 滤波器长度
    
    if L <= M && N > 0 % Allow L<=M if N=0 (empty input)
        % 在某些版本的重叠相加法中，L > M 是必需的，以确保 L-M+1 > 0
        % 对于本代码中的特定分段和累加逻辑 (valid_len = L-M+1 for num_segments and output_start)
        % 严格来说，如果L=M, valid_len=1. 如果L<M, valid_len<=0, num_segments 会非常大或负数.
        % 为了简单和鲁棒性，这里强制 L > M (除非输入为空)
        % 或者，修改分段逻辑以适应 L <= M 的情况（例如，标准重叠相加，输入分段长度L，输出累加间隔L）
        % 此处我们遵循用户原始代码的结构，其中 L-M+1 被用作计算分段数量和输出起始的有效长度。
        if L-M+1 <=0
           error('错误: 对于此overlapAdd实现, 分段长度 L (%d) 必须大于滤波器长度 M (%d)。', L, M);
        end
    end
    
    % 线性卷积结果长度
    P_out = N + M - 1;
    if N == 0 % Handle empty input signal
        y = zeros(P_out, 1); % Or zeros(0,1) if M=1, or based on conv behavior
        return;
    end

    % FFT 计算长度
    nfft = L + M - 1; 
    
    % 预计算 h 的 FFT (h 已经是列向量)
    H_fft = fft([h; zeros(nfft - M, 1)], nfft);
    
    % 每段有效输入/输出（非重叠部分）的长度，用于确定分段数和累加起始点
    % 这是用户原始代码的逻辑，我们将保留它。
    valid_len = L - M + 1;
    num_segments = ceil(N / valid_len);
        
    % 对输入信号x进行尾部补零，以确保可以取到所有分段
    % 需要补零的长度，使得 x_padded 的长度至少是 (num_segments-1)*valid_len + L
    required_x_len = (num_segments-1)*valid_len + L;
    if N < required_x_len
        x_padded = [x; zeros(required_x_len - N, 1)];
    else
        x_padded = x;
    end
    
    % 初始化输出数组
    y = zeros(P_out, 1);
    
    % 对每个分段进行处理
    for i = 1:num_segments
        % 提取当前分段 (长度为L)
        start_idx = (i-1) * valid_len + 1;
        % x_segment 必须是 L 长。如果 x_padded 不够长，则末尾补零。
        % (上面的 x_padded 逻辑已确保 x_padded(start_idx : start_idx+L-1) 是有效的)
        x_segment = x_padded(start_idx : start_idx + L - 1);
        
        % 计算分段与滤波器的线性卷积（使用FFT）
        X_segment_fft = fft(x_segment, nfft); % x_segment是Lx1, fft自动补零到nfft
        Y_segment_conv_fft = X_segment_fft .* H_fft;
        y_conv_segment = real(ifft(Y_segment_conv_fft)); % 结果长度为 nfft
        
        % 获取输出位置并累加
        output_start_idx = (i-1) * valid_len + 1;
        output_end_idx_ideal = output_start_idx + nfft - 1;
        
        % 确保不会写入超过 y 的分配长度 P_out
        actual_output_end_idx = min(output_end_idx_ideal, P_out);
        length_to_add = actual_output_end_idx - output_start_idx + 1;
        
        if length_to_add > 0
            y(output_start_idx : actual_output_end_idx) = ...
                y(output_start_idx : actual_output_end_idx) + y_conv_segment(1:length_to_add);
        end
    end
    
    % 截取最终需要的长度 (理论上 y 已经接近 P_out，但为确保，可以截断)
    if length(y) > P_out
        y = y(1:P_out);
    end
end
