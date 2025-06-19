# 基于IMU和磁力计的EKF姿态融合
本项目实现了基于扩展卡尔曼滤波（EKF）的 IMU 与磁力计数据融合，用于姿态估计，快速提供初始航向角，应用于姿态和航向参考系统（AHRS）中。通过对惯性测量单元 (IMU) 的加速度计与陀螺仪数据，以及磁力计数据进行滤波融合，提高航向角（Yaw）估计的鲁棒性与精度。
## 1. 磁力计软硬铁校准
### *MagCalib.m*
#### 文件格式

**输入数据**：`mag.log`
示例中MAG的输出频率为**75Hz**
格式：*[t, mx, my, mz]*
t  -- 时间戳 (**s**)
*mx, my, mz* -- 磁力计读数, 单位: **μT**

**输出数据**：`mag_calibed.txt`（校准后的磁力计数据）
保留原磁力计的输出频率**75Hz**以及**时间戳**
格式：*[t, mx, my, mz]*
t  -- 时间戳 (**s**)
*mx, my, mz* -- 磁力计读数, 单位: **μT**

**软铁校准矩阵**：`M_inv`
**硬铁偏置**：`c`
补偿公式为：**calib = M_inv * (raw - c);  % 先减偏置，再乘软铁校准矩阵**
注意：raw为为补偿的原始磁力计数据，calib为补偿校准之后的磁力计数据，且已经归一化，后续使用EKF进行融合请保留这两个参数

在使用罗盘传感器之前，需要对其进行校准以消除两个主要误差。一个是失调误差，这原本是由传感器和电路的失调误差引起的。另一个是标度误差。这两种误差都容易受到周围磁环境的干扰。例如，如果有一个x轴向的外部磁场施加到传感器上，就会给出外部x轴失调误差。同时，x轴标度也将与y轴和z轴不同。

通常用于校准磁传感器的方法是在xy平面上转动传感器绕圈，然后抽取数据。一个地点的地磁场强度是一个常数值，因此绘制的数据应该是一个圆；然而，事实上，我们将看到一个椭圆形，这意味着我们需要移动椭圆并重新缩放到以零为中心的圆。

上述2D校准方法有一些缺点，并且需要用加速器来测量其倾斜度。我们使用3D球面拟合方法来校准罗盘传感器。首先，我们需要将传感器旋转到x-y-z空间中的每个方向，并在3D坐标中绘制其值。然后我们需要使用最小平方误差(MSE)方法将数据拟合为椭球面。

[![](https://www.analog.com/cn/_/media/images/analog-dialogue/en/volume-53/number-1/articles/strapdown-inertial-navigation-system-based-on-an-imu-and-a-geomagnetic-sensor/234033-fig-01.jpg?w=900&rev=3388c0271a60461f96d9cad30ea90f9c&sc_lang=zh)](https://www.analog.com/cn/_/media/images/analog-dialogue/en/volume-53/number-1/articles/strapdown-inertial-navigation-system-based-on-an-imu-and-a-geomagnetic-sensor/234033-fig-01.jpg?w=900&rev=3388c0271a60461f96d9cad30ea90f9c&sc_lang=zh)

为了校准传感器，需要拉伸或压缩拟合的椭球面并将其移至以零为中心的球面上。这里使用矩阵奇异值分解(SVD)方法来进行这种校准。校准后的球体如下图所示。

[![](https://www.analog.com/cn/_/media/images/analog-dialogue/en/volume-53/number-1/articles/strapdown-inertial-navigation-system-based-on-an-imu-and-a-geomagnetic-sensor/234033-fig-02.png?w=900&rev=b1df667850aa4684a6805cb02c14df91&sc_lang=zh)](https://www.analog.com/cn/_/media/images/analog-dialogue/en/volume-53/number-1/articles/strapdown-inertial-navigation-system-based-on-an-imu-and-a-geomagnetic-sensor/234033-fig-02.png?w=900&rev=b1df667850aa4684a6805cb02c14df91&sc_lang=zh)

校准后可以看到，测得的磁场强度（球半径）几乎恒定不变。

[![](https://www.analog.com/cn/_/media/images/analog-dialogue/en/volume-53/number-1/articles/strapdown-inertial-navigation-system-based-on-an-imu-and-a-geomagnetic-sensor/234033-fig-03.png?w=900&rev=776e31f850cb4ffd8824ad276f07d730&sc_lang=zh)](https://www.analog.com/cn/_/media/images/analog-dialogue/en/volume-53/number-1/articles/strapdown-inertial-navigation-system-based-on-an-imu-and-a-geomagnetic-sensor/234033-fig-03.png?w=900&rev=776e31f850cb4ffd8824ad276f07d730&sc_lang=zh)

## 2.IMU & MAG EKF融合
### *IMU_MAG_EKF.m*
#### 文件格式
**输入数据**：`imu.log`
示例中IMU的输出频率为**100Hz**
格式：*[t, gx, gy, gz，ax, ay, az]*
t  -- 时间戳 (**s**)
*ax, ay, az* -- 加速度计读数, 单位: **g**
*gx, gy, gz *-- 角速度读数, 单位: **rad/s**

**输入数据**：`mag_data.txt`
示例中MAG的输出频率为**75Hz**
格式：*[t, mx, my, mz]*
t  -- 时间戳 (**s**)
*mx, my, mz* -- 磁力计读数, 单位: **μT**

#### EKF 初始化
**坐标系定义**：坐标系均采用机体系 **前+X / 右+Y / 上+Z**

**磁场参考**：m_n = [0.03,11.98,-58.03]'
请使用实际测得的磁北磁力计数据作为磁场参考

**重力参考**：g_n = [0 0 9.80665].'
请使用当地的重力数据作为重力参考

**磁力计更新开关量**：use_mag_update = 1
如果不需要磁力计参与EKF滤波量测更新，可以置为**0**

**初始化航向角yaw**
使用前1秒的平均磁力计数据计算初始航向角

**EKF状态向量选取**：x = [q0, q1, q2, q3, bgx, bgy, bgz]'
选取姿态四元数**q**以及陀螺仪的三轴零偏不稳定性**bg**作为状态量进行EKF滤波处理
不断估计输出最佳的姿态

#### 状态传播
使用陀螺仪数据进行积分，对状态进行预测
#### 重力加速度量测更新
依靠重力加速度可以约束roll和pitch两个姿态角
#### 磁力计量测更新
依靠磁力计可以约束yaw航向角
只要有磁力计数据可用，就进行磁力计的量测更新
