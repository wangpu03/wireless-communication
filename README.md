<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

# wireless-communication
无线信道matlab仿真等文件
### Clake and Jakes channel simulation

### rake 接收机仿真
* rake 接收机

### rayleigh channel simulation
* rayleigh 信号以及jake接收机仿真
* test_signal 信号生成仿真
* test_signal_FFT 信号Fourier变换仿真
* test_power 信号功率谱变换仿真
* single_sin_power 单个频率的正弦信号的功率以及功率谱的计算
* multi_sin_power 多个频率的正弦信号的功率以及功率谱的计算

----
### 能量信号和功率信号

根据信号可以用能量式或者功率式表示可以分为能量信号和功率信号

#### 能量信号
如各类瞬变信号，在非电量测量中，常将被测信号转换为电压或者电流信号来处理。电压信号工作在单位电阻\\(R = 1\\)上瞬时功率为\\(P(t) = x^2(t)/R = x^2(t)\\)。瞬时功率对时间的积分即是该时间内的能量。通常不考虑量纲，直接把信号的平方及其时间的积分分别称为信号的功率和能量。当\\(x(t)\\)满足
$$\int_{-\infty}^{+\infty}x^2(t)\,dt < \infty.$$
则信号的能量有限，称为能量有限信号，简称能量信号。其满足能量有限条件，实际是满足绝对可积条件。
我们定义信号\\(f(t)\\)的能量（归一化处理）：有电压\\(f(t)\\)或者电流\\(f(t)\\)在1\\(\Omega\\)电阻上消耗的能量：
$$E = \int_{-\infty}^{+\infty}x^2(t)\,dt, (注释：E= u*i=u^2/R=u^2).$$

#### 功率信号
如各种周期信号、常值信号、阶跃信号等，若\\(x(t)\\)在区间\\((-\infty,+\infty)\\)的能量无限，不满足上面的可积条件，但在有限区间\\((-T/2,+T/2)\\)满足平均功率有限的条件
$$\lim_{T \to +\infty} \frac{1}{T} \int_{-T/2}^{T/2} x^2(t)\,dt < +\infty$$
我们定义信号\\(f(t)\\)的平均功率，为电压\\(f(t)\\)在1\\(\Omega\\)电阻上消耗的平均功率（简称功率）：
$$\lim_{T \to +\infty} \frac{1}{T} \int_{-T/2}^{T/2} x^2(t)\,dt$$


