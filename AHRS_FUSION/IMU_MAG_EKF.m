%% imu_mag_fusion_ekf.m
% EKF 融合IMU与磁力计
% -------------------------------------------------------------------------
% 输入文件格式：
%   imu_data.txt : [t, gx, gy, gz，ax, ay, az]
%       t  -- 时间戳 (s)
%       ax,ay,az -- 加速度计读数, 单位:g
%       gx,gy,gz -- 角速度读数, 单位:rad/s
%
%   mag_data.txt : [t, mx, my, mz]
%       t  -- 时间戳 (s)
%       mx,my,mz -- 磁力计读数, 单位:μT
%
% 坐标系均采用机体系 前+X / 右+Y / 上+Z
% -------------------------------------------------------------------------
% 输出：
%   roll/pitch/yaw 姿态变化曲线
%% ------------------------ 读入数据 -----------------------------------
imu     = readmatrix('imu4.log');  
mag_raw = readmatrix('mag4.log');

t_imu  = imu(:,1);   gyro_b = imu(:,2:4);   acc_b = imu(:,5:7);
t_mag  = mag_raw(:,1);  
mag_b = (M_inv * (mag_raw(:,2:4).' - c))';

N      = length(t_imu);
N_mag  = length(t_mag);

use_mag_update = 1;

%% ------------------------ EKF 初始化 ---------------------------------
g_n = [0 0 9.80665].';         %重力参考
m_n = [0.03,11.98,-58.03]';    %磁场参考(朝磁北)
m_n = m_n / norm(m_n);

idx_ref = t_mag <= (t_mag(1) +1);
m_s = mean(mag_b(idx_ref,:),1).';  
m_s = m_s / norm(m_s);
yaw_rad = atan2(m_s(1), m_s(2)); 
quat = [cos(yaw_rad / 2), 0, 0, sin(yaw_rad / 2)];
quat = quat / norm(quat);

% 四元数初始全部 1 0 0 0
q_arr  = quaternion(ones(N,1), zeros(N,1), zeros(N,1), zeros(N,1));
x  = zeros(7,1);

x(1:4) = quat.';               % 四元数(初值 前1秒钟平均yaw)
x(5:7) = [0;0;0];              % 陀螺仪漂移 (初值 0)

P  = blkdiag(eye(4)*1e-3, eye(3)*(0.5*pi/180)^2);   % 协方差

%% 噪声设置
sigma_gyro = 0.02*pi/180;      % 角速度噪声 (rad/s)
sigma_bg   = 0.001*pi/180;     % 漂移随机游走 (rad/s)
Q_gyro = (sigma_gyro)^2 * eye(3);
Q_bg   = (sigma_bg)^2   * eye(3);

R_acc  = (0.03*9.80665)^2 * eye(3);
R_mag  = (0.06)^2        * eye(3);

%% ------------------------ 主循环 -------------------------------------
skew = @(v)[  0   -v(3)  v(2);
              v(3)   0  -v(1);
             -v(2)  v(1)   0 ];

idxMag = 1;
tolMag = 0.002;  

eps = 1e-7;     

for k = 2:N
    dt = t_imu(k) - t_imu(k-1);  if dt<=0, dt=1/100; end
    
    omega_m = gyro_b(k,:).';          % 测量角速度
    b_g     = x(5:7);                 % 漂移估计
    omega   = omega_m - b_g;          % 去漂后的真实角速度
    
    Omega = [ 0      -omega.' ;
              omega   -skew(omega) ]; % 4×4
    
    % 四元数积分 (一阶)
    x_pred           = x;
    x_pred(1:4)      = x(1:4) + 0.5*Omega*x(1:4)*dt;
    x_pred(1:4)      = x_pred(1:4) / norm(x_pred(1:4));
    x_pred(5:7)      = b_g;           % 漂移为常值
    
    % 状态转移雅可比 F (7×7)
    F      = eye(7);
    F(1:4,1:4) = eye(4) + 0.5*Omega*dt;
    Gq     = 0.5*dt * [ -x(2:4).';
                        x(1)*eye(3) + skew(x(2:4)) ];  % 4×3
    F(1:4,5:7) = -Gq;               % ∂q/∂b_g
    
    % 噪声耦合矩阵 G (7×6)
    G = zeros(7,6);
    G(1:4,1:3) = Gq;      % gyro 噪声
    G(5:7,4:6) = eye(3);  % 漂移随机游走
    
    Qd = blkdiag(Q_gyro, Q_bg);      % 6×6
    
    P = F*P*F.' + G*Qd*G.';          % 协方差传播
    
    %% ========= 加计更新 =========
    % 预测重力
    q_pred   = quaternion(x_pred(1:4).');
    R_bn     = rotmat(q_pred,'frame');
    g_pred   = R_bn * g_n;
    
    % 雅可比 H_a (3×7)
    H_q = zeros(3,4);
    for i = 1:4
        dq  = zeros(4,1); dq(i)=eps;
        q_e = x_pred(1:4)+dq;  q_e = q_e / norm(q_e);
        R_e = rotmat(quaternion(q_e.'),'frame');
        H_q(:,i) = (R_e*g_n - g_pred)/eps;
    end
    H_a = [H_q, zeros(3,3)];
    
    z_acc = acc_b(k,:).' * 9.80665;           % g → m/s^2
    y_a   = z_acc - g_pred;
    S_a   = H_a*P*H_a.' + R_acc;
    K_a   = P*H_a.'/S_a;
    x_upd = x_pred + K_a*y_a;
    P     = (eye(7)-K_a*H_a)*P;
    
    %% ========= 磁计更新 =========
    if(use_mag_update == 0)
        x_upd(1:4) = x_upd(1:4) / norm(x_upd(1:4));
        x = x_upd;
        q_arr(k) = quaternion(x(1:4).');
    else
        while idxMag <= N_mag && (t_mag(idxMag) < t_imu(k)+tolMag)
            if abs(t_mag(idxMag)-t_imu(k)) <= tolMag
                z_mag = mag_b(idxMag,:).';
                q_upd = quaternion(x_upd(1:4).');
                R_bn  = rotmat(q_upd,'frame');
                m_pred = R_bn * m_n;

                % 雅可比 H_m (3×7)
                H_q = zeros(3,4);
                for i = 1:4
                    dq  = zeros(4,1); dq(i)=eps;
                    q_e = x_upd(1:4)+dq;  q_e = q_e / norm(q_e);
                    R_e = rotmat(quaternion(q_e.'),'frame');
                    H_q(:,i) = (R_e*m_n - m_pred)/eps;
                end
                H_m = [H_q, zeros(3,3)];

                y_m = z_mag - m_pred;
                S_m = H_m*P*H_m.' + R_mag;
                K_m = P*H_m.'/S_m;
                x_upd = x_upd + K_m*y_m;
                P     = (eye(7)-K_m*H_m)*P;
            
            end
            idxMag = idxMag + 1;
        end
        % 归一化 & 保存
        x_upd(1:4) = x_upd(1:4) / norm(x_upd(1:4));
        x = x_upd;
        q_arr(k) = quaternion(x(1:4).');
    end
end

%% ------------------------ 作图 ---------------------------
eul = eulerd(q_arr(2:k),'ZYX','frame');  
figure;
subplot(3,1,1); plot(t_imu(2:k),eul(:,3)); ylabel('Roll (°)');
subplot(3,1,2); plot(t_imu(2:k),eul(:,2)); ylabel('Pitch (°)');
subplot(3,1,3); plot(t_imu(2:k),eul(:,1)); ylabel('Yaw (°)');
xlabel('Time (s)'); sgtitle('EKF 姿态估计');

q = compact(q_arr); 
fig = figure('Name','Attitude 3‑D Animation');
ax  = axes('Parent',fig);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); view(ax,3);
xlabel(ax,'X_N'); ylabel(ax,'Y_N'); zlabel(ax,'Z_N');
axis(ax,[-1.2 1.2 -1.2 1.2 -1.2 1.2]);

% 画一个透明立方体代表机体
[V,F] = createCube(0.2);         % 立方体边长 0.2
patch('Vertices',V,'Faces',F,'FaceColor',[0.8 0.8 1],...
      'FaceAlpha',0.3,'EdgeColor','k','Parent',ax);

% 三根机体坐标轴 (初始全在 +方向)
bx = quiver3(ax,0,0,0,1,0,0,'r','LineWidth',2,'MaxHeadSize',0.5);
by = quiver3(ax,0,0,0,0,1,0,'g','LineWidth',2,'MaxHeadSize',0.5);
bz = quiver3(ax,0,0,0,0,0,1,'b','LineWidth',2,'MaxHeadSize',0.5);

% 时间文本
tt = text(ax,0.05,0.05,1.1,'','Units','normalized');

%% 3. 动画循环
skip = max(1,floor(N/500));   % 若样本过多，跳帧避免过慢
for k = 1:skip:N
    R = quat2dcm(q(k,:)); % 3×3, Nav→Body
    % 机体三轴在导航下
    ex = R(:,1); ey = R(:,2); ez = R(:,3);
    set(bx,'UData',ex(1),'VData',ex(2),'WData',ex(3));
    set(by,'UData',ey(1),'VData',ey(2),'WData',ey(3));
    set(bz,'UData',ez(1),'VData',ez(2),'WData',ez(3));
    % 更新立方体姿态
    Vb = (R * V.').';           % 把顶点由 Body→Nav
    set(findobj(ax,'Type','patch'),'Vertices',Vb);

    set(tt,'String',sprintf('t = %.2f s',t_imu(k)));
    drawnow;
    if ~ishandle(fig); break; end % 关闭窗口则退出
end

function [V,F] = createCube(a)
    % a: 边长 (默认 1)
    if nargin<1; a = 1; end
    A = a/2;
    V = [ -A -A -A;
           +A -A -A;
           +A +A -A;
           -A +A -A;
           -A -A +A;
           +A -A +A;
           +A +A +A;
           -A +A +A ];
    F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
end
