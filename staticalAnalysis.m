
% 使用滑动窗口截取分析信号片段特征
% 先运行BWN_1
signal = ECG_lp_filtered; 

window_size = 150;  % 滑动窗口长度
n_windows = length(signal) - window_size + 1;

% 储峰度、偏度、标准差
kurtosis_vals = zeros(1, n_windows);
skewness_vals = zeros(1, n_windows);
std_vals = zeros(1, n_windows);

for i = 1:n_windows
    window = signal(i:i+window_size-1); 
    kurtosis_vals(i) = kurtosis(window); 
    skewness_vals(i) = skewness(window); 
    std_vals(i) = std(window); 
end

subplot(4,1,1);
plot(signal);
title('Original Signal');

subplot(4,1,2);
plot(kurtosis_vals);
title('Kurtosis ');

subplot(4,1,3);
plot(skewness_vals);
title('Skewness');

subplot(4,1,4);
plot(std_vals);
title('Standard Deviation');



std_threshold = 50;  % 标准差阈值
kurtosis_threshold = 10;  

% 标准差超过阈值的索引
outlier_indices_std = find(std_vals > std_threshold);
% outlier_indices_kurtosis = find(kurtosis_vals > kurtosis_threshold);


% outlier_indices = unique([outlier_indices_std, outlier_indices_kurtosis]);
outlier_indices = outlier_indices_std;

% 起始和结束索引
start_index = min(outlier_indices);
end_index = max(outlier_indices);

disp(['Detected disturbance between indices ', num2str(start_index), ' and ', num2str(end_index)]);

% 对区间信号进行移动平均滤波
smoothed_signal = signal; 
window_size = 300;  % 移动平均滤波窗口大小

smoothed_signal(start_index:end_index) = movmean(signal(start_index:end_index), window_size);

figure;
subplot(2,1,1);
plot(signal);
title('Original Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2);
plot(smoothed_signal);
title('Smoothed Signal');
xlabel('Time');
ylabel('Amplitude');
