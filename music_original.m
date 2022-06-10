clear;clc;
% ���������MUSIC��MVDR��root-MUSIC

%% ��������
N = 12;                              % ��Ԫ����
M = 7;                              % ��Դ����
theta = [-55 -40 -20 0 20 40 55];                 % �ź�����(��)
snr = 3;                            % �����
K = 1024;                           % ������

%% ����M����������������ź���λ��
wavelen = 1;                                % �����źŲ���Ϊ1
dd = wavelen / 2;                           % ��Ԫ���d = wavelen/2
d = 0:dd:(N-1)*dd;                          % ������������
A = exp(1i.*2*pi*d.'*sind(theta)/wavelen);  % �����������Σ����ź�����

%% ����M������źţ�����1024
S = randn(M, K) + 1i .* randn(M, K);        % ����������ź�

X = A * S;                                  % ���ź�������з���
X1 = awgn(X, snr, 'measured');              % ��������������С����ֵ

%% �����źŵ�Э�������
Rxx = X1 * X1' / K;                         % ����Э�������
[EV, D] = eig(Rxx);                         % �õ�������EV + ����ֵD  �°汾matlab�Ѿ���С�����������
figure(7)
bar3(D);
title('����ֵ��������')
% Rxx = [Rxx(:,3:end) Rxx(:,1:2)];

%% ѭ��������������������ʱ��
tic
idx = 1;
scale = -90 : 0.1 : 90;
[SP, SP_inv] = deal(zeros(length(scale), 1));
mvdr = zeros(length(scale), 1);    
En = EV(:, 1:end-M);                                    % С����ֵ����������
invRxx = Rxx^(-1);
for angle_degree = scale
    a = exp(1j .* 2*pi*d*sind(angle_degree)/wavelen).'; % �����źŵ���ʸ�����ù���ת��ȫ���Ӹ���                           
    SP(idx) = (a'*En)*(En'*a);                          % ����ǰ�潲���������жϽ��
    SP_inv(idx) = 1/abs((a'*En)*(En'*a));               % �õ�����תһ�� ��ɷ�ֵ
    
    mvdr(idx) = 1/abs(a'*(invRxx)*a);
    
    idx = idx + 1;
end
% SP_db = db(SP);                                         % ת��ΪdB
% SP_inv_db = db(abs(SP_inv));                            % ת��ΪdB

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
z = roots(co);                  % ����ʽ���

rx = z.';
[as, ad] = sort(abs(abs(rx) - 1));
doa_rmusic = asin(sort(-angle(rx(ad([1:2:2*M]))) / pi)) * 180/pi
toc

%% ploting
figure(8)
% subplot(211);
% plot(scale, SP_db);
% xlabel('�����/(degree)'); ylabel('�ռ���/(dB)');
% grid on;
% title('������ʾĿ���ź�����')
% 
% subplot(212);
plot(scale, SP_invdb, 'b');
xlabel('�����/(degree)'); ylabel('�ռ���/(dB)');
grid on;
% title('��MUSIC�ױ�ʾĿ���ź�����')
hold on;
plot(scale, mvdr_db, 'r');
xlabel('�����/(degree)'); ylabel('�ռ���/(dB)');
grid on;
% title('��Capon(MVDR)�ױ�ʾĿ���ź�����')
legend('MUSIC', 'MVDR');
title('Ŀ���ź�����DOA');

% ���ߣ�qwe14789cn
% ���ӣ�https://zhuanlan.zhihu.com/p/362008506
% ��Դ��֪��
% ����Ȩ���������С���ҵת������ϵ���߻����Ȩ������ҵת����ע��������

% Ϊ��d<lambda/2
% https://www.zhihu.com/question/265573245/answer/1945780095

