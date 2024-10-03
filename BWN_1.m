clear all;
load("ECG_database.mat");

%% RLS
fc = 0.5;  
[b_lp, a_lp] = butter(2, fc/(FS/2), 'low');
reference_signal = filtfilt(b_lp, a_lp, BWN_data);

reference_signal = reference_signal(:);

M = 1;             % 滤波器阶数
delta = 0.1;       % 初始化P的小正数
lambda = 0.99;     % 遗忘因子
P = (1/delta)*eye(M);
w = zeros(M,1);

y = zeros(LENGTH,1);
e = zeros(LENGTH,1);

for n = 1:LENGTH
    x_n = reference_signal(n);

    d_n = BWN_data(n);

    y_n = w' * x_n;

    e_n = d_n - y_n;

    k_n = (P * x_n) / (lambda + x_n' * P * x_n);

    w = w + k_n * e_n;

    P = (1/lambda)*(P - k_n * x_n' * P);

    y(n) = y_n;
    e(n) = e_n;
end

%% 高通低通
fc_high = 1; 
[b_hp, a_hp] = butter(4, fc_high/(FS/2), 'high');  
ECG_hp_filtered = filtfilt(b_hp, a_hp, e); 

fc_low = 40;  
[b_lp, a_lp] = butter(4, fc_low/(FS/2), 'low'); 
ECG_lp_filtered = filtfilt(b_lp, a_lp, ECG_hp_filtered);  

%% 
figure;
subplot(3,1,1);
plot(BWN_data);
title('BWNECG信号');

subplot(3,1,2);
plot(e);
title('RLS滤波');

subplot(3,1,3);
plot(ECG_lp_filtered);
title('高通低通滤波后');

