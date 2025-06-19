% 磁力计数据椭球拟合与可视化
clear; clc; close all;

% 文件读取 - 健壮的数据处理方法
filename = 'mag4.log';
fprintf('正在读取文件: %s\n', filename);

try
    % 尝试多种读取方法
    data = readmatrix(filename);
    
    % 如果出现NaN，尝试跳过标题行
    if any(isnan(data(:)))
        fprintf('检测到NaN值，尝试跳过标题行...\n');
        opts = detectImportOptions(filename);
        opts.DataLines = [2, Inf];
        data = readmatrix(filename, opts);
    end
    
    % 如果仍有问题，使用文本扫描
    if any(isnan(data(:))) || isempty(data)
        fprintf('使用文本扫描读取数据...\n');
        fid = fopen(filename, 'r');
        rawData = textscan(fid, '%f %f %f', 'HeaderLines', 1, 'Delimiter', ',; \t', 'MultipleDelimsAsOne', true);
        fclose(fid);
        data = [rawData{1}, rawData{2}, rawData{3}];
    end
    
    
    % 提取XYZ分量
    x = data(:, 2);
    y = data(:, 3);
    z = data(:, 4);
    
    % 移除任何NaN值
    validIdx = ~isnan(x) & ~isnan(y) & ~isnan(z);
    x = x(validIdx);
    y = y(validIdx);
    z = z(validIdx);
    
    fprintf('成功读取 %d 个有效数据点\n', length(x));
    
    % 检查是否有足够的数据点
    if length(x) < 10
        error('数据点不足，至少需要10个点进行椭球拟合');
    end
    
catch ME
    fprintf('文件读取错误: %s\n', ME.message);
    return;
end

% 创建图形
figure('Color', 'white', 'Name', '磁力计数据椭球拟合', 'Position', [100, 100, 1000, 800]);

% 绘制原始数据点
scatter3(x, y, z, 50, 'filled', ...
         'MarkerFaceColor', [0.2, 0.6, 0.9], ...
         'MarkerEdgeColor', 'k', ...
         'MarkerFaceAlpha', 0.7);
hold on;
grid on;

% =============================================
% 最小二乘法椭球拟合（修正版）
% =============================================
fprintf('正在进行最小二乘椭球拟合...\n');

% 椭球一般方程: 
% A*x^2 + B*y^2 + C*z^2 + D*x*y + E*x*z + F*y*z + G*x + H*y + I*z + J = 0

% 构建设计矩阵
X = [x.^2, y.^2, z.^2, x.*y, x.*z, y.*z, x, y, z, ones(size(x))];

% 约束条件：A + B + C = 3（防止零解）
C = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0];  % 约束向量

% 构建增广矩阵
A = X' * X;
M = [2*A, C; C', 0];

% 右侧向量 (零向量 + 约束值)
b = [zeros(10,1); 3];

% 求解线性系统
params = M \ b;
coeffs = params(1:10);  % 提取椭球系数

% 提取系数
A = coeffs(1); B = coeffs(2); C = coeffs(3);
D = coeffs(4); E = coeffs(5); F = coeffs(6);
G = coeffs(7); H = coeffs(8); I = coeffs(9);
J = coeffs(10);

fprintf('椭球拟合完成\n系数: A=%.4f, B=%.4f, C=%.4f, D=%.4f, E=%.4f, F=%.4f, G=%.4f, H=%.4f, I=%.4f, J=%.4f\n', ...
        A, B, C, D, E, F, G, H, I, J);

% =============================================
% 将一般椭球方程转换为标准形式
% =============================================

% 构造二次型矩阵
Q = [A, D/2, E/2;
     D/2, B, F/2;
     E/2, F/2, C];
 
% 构造线性项向量
T = [G; H; I];

% 计算椭球中心
center = -Q \ T / 2;
fprintf('椭球中心: (%.4f, %.4f, %.4f)\n', center(1), center(2), center(3));

% 计算平移后的常数项
s = center' * Q * center + T' * center + J;

% 计算特征值和特征向量
[R, S] = eig(Q / (-s));  % Q/s 的特征分解

% 计算半轴长度
radii = 1 ./ sqrt(diag(S));
fprintf('椭球半轴长度: a=%.4f, b=%.4f, c=%.4f\n', radii(1), radii(2), radii(3));

% 生成椭球表面
[ex, ey, ez] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 50);

% 旋转椭球到主方向
for i = 1:size(ex,1)
    for j = 1:size(ex,2)
        vec = [ex(i,j); ey(i,j); ez(i,j)];
        vec_rot = R * vec;
        ex(i,j) = vec_rot(1) + center(1);
        ey(i,j) = vec_rot(2) + center(2);
        ez(i,j) = vec_rot(3) + center(3);
    end
end

% 绘制拟合椭球
surf(ex, ey, ez, 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
     'FaceColor', [0.9, 0.4, 0.2], 'FaceLighting', 'gouraud');
 
% =============================================
% 可视化设置
% =============================================

% 设置坐标轴和标签
box on;
axis equal;
xlabel('X轴场强 (μT)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Y轴场强 (μT)', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('Z轴场强 (μT)', 'FontSize', 14, 'FontWeight', 'bold');
title('椭球拟合', 'FontSize', 16, 'FontWeight', 'bold');

% 添加坐标轴指示线
line_x = line([min(x), max(x)], [0, 0], [0, 0], 'Color', 'r', 'LineWidth', 2.5);
line_y = line([0, 0], [min(y), max(y)], [0, 0], 'Color', [0, 0.7, 0], 'LineWidth', 2.5);
line_z = line([0, 0], [0, 0], [min(z), max(z)], 'Color', 'b', 'LineWidth', 2.5);

% 添加原点标记
plot3(0, 0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');

% 添加椭球中心标记
plot3(center(1), center(2), center(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% 添加光照效果
light('Position', [1, 1, 1], 'Style', 'infinite', 'Color', [1, 1, 0.8]);
light('Position', [-1, -1, -1], 'Style', 'infinite', 'Color', [0.8, 0.8, 1]);
lighting gouraud;
material dull;

% 添加图例
legend_items = {'原始数据点', '拟合椭球', 'X轴', 'Y轴', 'Z轴', '坐标原点', '椭球中心'};
legend(legend_items, 'Location', 'best');

% 计算拟合误差
F2 = A*x.^2 + B*y.^2 + C*z.^2 + D*x.*y + E*x.*z + F*y.*z + G*x + H*y + I*z + J;
mse = mean(F2.^2);
fprintf('均方误差(MSE): %.6f\n', mse);

% 添加数据统计信息
stats = sprintf(['NUM: %d\n'...
                'X: [%.2f, %.2f] μT\n'...
                'Y: [%.2f, %.2f] μT\n'...
                'Z: [%.2f, %.2f] μT\n'...
                'Fitting error(MSE): %.4f μT\n'...
                'Ellipsoid center: (%.2f, %.2f, %.2f)\n'...
                'Axis length: %.2f, %.2f, %.2f'], ...
                length(x), min(x), max(x), min(y), max(y), min(z), max(z), ...
                mse, ...
                center(1), center(2), center(3), ...
                radii(1), radii(2), radii(3));
            
annotation('textbox', [0.75, 0.15, 0.2, 0.25], 'String', stats, ...
           'FitBoxToText', 'on', 'BackgroundColor', [1, 1, 1, 0.7], ...
           'EdgeColor', 'black', 'FontSize', 11, 'FontName', 'Consolas');

% 设置视角
view(135, 30);
rotate3d on;

% 添加坐标轴箭头
arrow_start = [min(x) min(y) min(z)];
arrow_end = [max(x) max(y) max(z)];

% X轴箭头
quiver3(0, 0, 0, max(x)*1.1, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% Y轴箭头
quiver3(0, 0, 0, 0, max(y)*1.1, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% Z轴箭头
quiver3(0, 0, 0, 0, 0, max(z)*1.1, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);

params = [A,B ,C ,D, E, F ,G, H, I, J];
a = [ params(1)  params(4)/2  params(5)/2;
      params(4)/2  params(2)  params(6)/2;
      params(5)/2  params(6)/2  params(3) ];
b = params(7:9).';
d = params(10);

% 1) 计算硬铁偏置
c = -0.5 * (a \ b);     % 等价于 inv(A)*b

% 2) 归一化
k  = c.'*a*c - d;
A_ = a / k;

% 3) SVD 分解
[U,S,~] = svd(A_);      % A_ = U*S*U'
s = diag(S);            % 对角元素 σ_i

% 4) 构造校准矩阵 (这里取乘以 M^(-1) 的形式更常用)
M_inv = U * diag(sqrt(s)) * U';    % 软铁补偿

%% 对一批原始数据 [x y z] 做校准
raw   = [x y z].';          % 3×N
calib = M_inv * (raw - c);  % 先减偏置，再乘软铁校准矩阵

%% 可视化验证

% ---------- 2. 计算半径 ----------
radius_pre  = sqrt(sum(data(:,2:4).^2 , 2));   % Nx1
R0          = mean(radius_pre);       
caca = calib';
caca = caca * R0;                   

figure; scatter3(caca(:,1), caca(:,2), caca(:,3), 5, 'filled',...
         'MarkerFaceColor', [0.2, 0.6, 0.9], ...
         'MarkerEdgeColor', 'k', ...
         'MarkerFaceAlpha', 0.7);
axis equal; grid on;

% 设置坐标轴和标签
box on;
axis equal;
xlabel('X轴场强 (μT)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Y轴场强 (μT)', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('Z轴场强 (μT)', 'FontSize', 14, 'FontWeight', 'bold');

% 添加坐标轴指示线
line_x = line([min(caca(:,1)), max(caca(:,1))], [0, 0], [0, 0], 'Color', 'r', 'LineWidth', 2.5);
line_y = line([0, 0], [min(caca(:,2)), max(caca(:,2))], [0, 0], 'Color', [0, 0.7, 0], 'LineWidth', 2.5);
line_z = line([0, 0], [0, 0], [min(caca(:,3)), max(caca(:,3))], 'Color', 'b', 'LineWidth', 2.5);

radius_post = sqrt(sum(caca.^2, 2));   % Nx1

% ---------- 3. 绘图 ----------
pts = (1:numel(radius_pre)).';           % 横轴 = 点编号

figure('Color','w'); hold on;
plot(pts, radius_pre,  'Color',[0.2 0.6 1 0.7], 'LineWidth',1);    % 校准前
plot(pts, radius_post, 'Color',[0.6 0.1 0.1], 'LineWidth',1.5);    % 校准后

title('\bf Magnetic Field Strength (Digitalized)');
xlabel('Points'); ylabel('Radius (uT)');
legend({'Pre-calibration','Post-calibration'},'Location','southwest');
grid on; box on; axis tight;
set(gca,'FontName','Helvetica','LineWidth',1.2,'YMinorGrid','on','XMinorTick','on');
fprintf('可视化完成\n');


out = [ data(:,1) , caca ];   % Nx4   [timestamp, Bx, By, Bz]
writematrix(out, 'mag_calibed.txt', 'Delimiter','tab');  
fprintf('文件写入完成\n');


