
clear all;
load("ECG_database.mat");

% 将PLI噪声添加到ECG数据
PLI_data = PLI_data + 20 * mains_signal;

fs = FS;

N = LENGTH;

t = (0:N-1) / fs;

%% RLS

% 50 Hz的正弦和余弦分量
x_ref = [sin(2*pi*50*t); cos(2*pi*50*t)]';  % N x 2矩阵

M = 2;               % 滤波器阶数（对应于正弦和余弦分量）
delta = 0.1;         % 初始化P的小正数
lambda = 0.99;       % 遗忘因子（接近1）
P = (1/delta)*eye(M);% 相关矩阵的逆
w = zeros(M,1);      % 初始滤波器权重（列向量）

y_rls = zeros(N,1);      % 滤波器输出（估计的PLI噪声）
e_rls = zeros(N,1);      % 误差信号（滤波后的ECG）

% RLS算法
for n = 1:N
    % 输入向量（第n时刻的参考信号）
    x_n = x_ref(n,:)';   % M x 1列向量

    % 期望信号（第n时刻的含噪声ECG）
    d_n = PLI_data(n);

    % 滤波器输出（估计的PLI噪声）
    y_n = w' * x_n;

    % 误差信号（期望信号与估计噪声之差）
    e_n = d_n - y_n;

    % 增益向量
    k_n = (P * x_n) / (lambda + x_n' * P * x_n);

    % 更新滤波器权重
    w = w + k_n * e_n;

    % 更新相关矩阵的逆
    P = (1/lambda)*(P - k_n * x_n' * P);

    % 存储输出
    y_rls(n) = y_n;
    e_rls(n) = e_n;
end

%% 使用Butterworth带阻滤波器去除PLI噪声

% 定义带阻滤波器的停止带频率
f_stop1 = 48;   % 下限频率（Hz）
f_stop2 = 52;   % 上限频率（Hz）

% 计算归一化频率（0到1之间）
fn = fs / 2;                   % 奈奎斯特频率
w_stop1 = f_stop1 / fn;        % 归一化下限频率
w_stop2 = f_stop2 / fn;        % 归一化上限频率

% 设计Butterworth带阻滤波器
filter_order = 5; % 滤波器阶数（可根据需要调整）
[b, a] = butter(filter_order, [w_stop1 w_stop2], 'stop');

% 使用filtfilt函数对数据进行零相位滤波
filtered_data = filtfilt(b, a, PLI_data);

%% 绘制时域信号的比较

figure;
plot(t, Data1, 'k', 'DisplayName', '无PLI噪声原始信号'); hold on;
plot(t, e_rls, 'r', 'DisplayName', 'RLS滤波后信号');
plot(t, filtered_data, 'b', 'DisplayName', '带阻滤波器滤波后信号');
legend('show');
xlabel('时间 (s)');
ylabel('幅度');
title('时域信号的比较');
grid on;

%% 计算原始和滤波后信号的FFT

Y_orig = fft(PLI_data);          % 原始信号的FFT
Y_rls = fft(e_rls);              % RLS滤波后信号的FFT
Y_butter = fft(filtered_data);   % 带阻滤波器滤波后信号的FFT

% 生成用于绘图的频率向量（直到奈奎斯特频率）
f = (0:N/2-1) * (fs / N);        % 频率向量，范围为0到fs/2

% 计算幅度谱（正频率部分）
Y_orig_mag = abs(Y_orig(1:N/2));   % 原始信号的幅度谱
Y_rls_mag = abs(Y_rls(1:N/2));     % RLS滤波后信号的幅度谱
Y_butter_mag = abs(Y_butter(1:N/2)); % 带阻滤波器滤波后信号的幅度谱

%% 绘制频域信号的比较

figure;
plot(f, Y_orig_mag, 'k', 'DisplayName', '含PLI噪声频谱');hold on;
plot(f, Y_rls_mag, 'r', 'DisplayName', 'RLS滤波后频谱');  
plot(f, Y_butter_mag, 'b', 'DisplayName', '带阻滤波器滤波后频谱');
legend('show');
xlabel('频率 (Hz)');
ylabel('幅度');
title('频域信号的比较');
grid on;
