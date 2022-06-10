clear;clc;
% 本代码包含MUSIC、MVDR、root-MUSIC

%% 参数设置
N = 12;                              % 阵元个数
M = 7;                              % 信源个数
theta = [-55 -40 -20 0 20 40 55];                 % 信号来向(°)
snr = 3;                            % 信噪比
K = 1024;                           % 快拍数

%% 生成M个来向的阵列流形信号相位差
wavelen = 1;                                % 来波信号波长为1
dd = wavelen / 2;                           % 阵元间距d = wavelen/2
d = 0:dd:(N-1)*dd;                          % 构建阵列坐标
A = exp(1i.*2*pi*d.'*sind(theta)/wavelen);  % 构建阵列流形，即信号来向

%% 发射M组随机信号，长度1024
S = randn(M, K) + 1i .* randn(M, K);        % 构建不相关信号

X = A * S;                                  % 对信号来向进行仿真
X1 = awgn(X, snr, 'measured');              % 加入噪声，产生小特征值

%% 构建信号的协方差矩阵
Rxx = X1 * X1' / K;                         % 构建协方差矩阵
[EV, D] = eig(Rxx);                         % 拿到特向量EV + 特征值D  新版本matlab已经从小到大排序好了
figure(7)
bar3(D);
title('特征值矩阵排列')
% Rxx = [Rxx(:,3:end) Rxx(:,1:2)];

%% 循环搜索特征向量正交的时候
tic
idx = 1;
scale = -90 : 0.1 : 90;
[SP, SP_inv] = deal(zeros(length(scale), 1));
mvdr = zeros(length(scale), 1);    
En = EV(:, 1:end-M);                                    % 小特征值的特征向量
invRxx = Rxx^(-1);
for angle_degree = scale
    a = exp(1j .* 2*pi*d*sind(angle_degree)/wavelen).'; % 构建信号导向矢量，用共轭转至全部加负号                           
    SP(idx) = (a'*En)*(En'*a);                          % 利用前面讲的正交来判断结果
    SP_inv(idx) = 1/abs((a'*En)*(En'*a));               % 用倒数翻转一下 变成峰值
    
    mvdr(idx) = 1/abs(a'*(invRxx)*a);
    
    idx = idx + 1;
end
% SP_db = db(SP);                                         % 转换为dB
% SP_inv_db = db(abs(SP_inv));                            % 转换为dB

SP_inv_max = max(SP_inv);
SP_invdb = 10 * log10(SP_inv / SP_inv_max);

mvdr_max = max(mvdr);
mvdr_db = 10 * log10(mvdr / mvdr_max);
toc

%% root-music
tic
GG = En * En';
co = zeros(2*N-1, 1);        
for m = 1:N
    co(m:m+N-1) = co(m:m+N-1) + GG(N:-1:1,m); 
end
z = roots(co);                  % 多项式求根

rx = z.';
[as, ad] = sort(abs(abs(rx) - 1));
doa_rmusic = asin(sort(-angle(rx(ad([1:2:2*M]))) / pi)) * 180/pi
toc

%% ploting
figure(8)
% subplot(211);
% plot(scale, SP_db);
% xlabel('入射角/(degree)'); ylabel('空间谱/(dB)');
% grid on;
% title('正交表示目标信号来向')
% 
% subplot(212);
plot(scale, SP_invdb, 'b');
xlabel('入射角/(degree)'); ylabel('空间谱/(dB)');
grid on;
% title('用MUSIC谱表示目标信号来向')
hold on;
plot(scale, mvdr_db, 'r');
xlabel('入射角/(degree)'); ylabel('空间谱/(dB)');
grid on;
% title('用Capon(MVDR)谱表示目标信号来向')
legend('MUSIC', 'MVDR');
title('目标信号来向DOA');

% 作者：qwe14789cn
% 链接：https://zhuanlan.zhihu.com/p/362008506
% 来源：知乎
% 著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。

% 为何d<lambda/2
% https://www.zhihu.com/question/265573245/answer/1945780095

