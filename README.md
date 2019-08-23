# Gammatone-filters
Python implementation of all-pole Gammatone filter
## Algorithms
$g(t)$ 可分解为两部分的乘积，即
$$
g(t)=a \times r(t) \times s(t)
$$
其中
$$
\begin{align}
r(t)&=t^{n-1}e^{-2\pi bt}\\
s(t)&=cos(2\pi f_c t+\phi)
\end{align}
$$

时域相乘==频域卷积，即：
$$
G(f)=a\times R(f)*S(f)
$$

可以分别计算 $R(f)$ 和 $S(f)$ ，即：
$$
\begin{equation}
\begin{aligned}
R(f)=FT(t^{n-1}e^{-2\pi b t})
&=\frac{1}{(j2\pi)^{n-1}}\frac{\partial^{n-1} FT(e^{-2\pi bt})}{\partial f^{n-1}}\\
&=\frac{1}{(j2\pi)^{n-1}}\frac{\partial^{n-1}\frac{1}{2\pi b+j2\pi f}}{\partial f^{n-1}}\\
&=\frac{1}{(j2\pi)^{n-1}}\frac{(j)^{n-1}(n-1)!}{2\pi}\frac{1}{(b+jf)^n}\\
&=\frac{(n-1)!}{(2\pi b)^n}\frac{1}{(1+jf/b)^n}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
S(f)=FT\left(cos(2\pi f_ct+\phi)\right)
&=e^{j\phi}\delta(f-f_c)+e^{-j\phi}\delta(f+f_c)
\end{aligned}
\end{equation}
$$

$R(f)$ 是低通滤波器，由n个一阶低通滤波器及联得到；$S(f)$ 的功能则是频率搬移，将低通滤波器转换为带通滤波器。
在设计滤波器的时候，可以反过来，首先对输入信号降频 $f_c$ ，然后在应用低通滤波器 $R(f)$ ，问题就变的简单了。
1. 降频

  $$
  \begin{equation}
  \begin{aligned}
  x'(t)=IFT(X(f-f_c))&=IFT\left(\int{x(t)e^{-j2\pi(f-f_c)t}dt}\right)\\
  &=IFT\left(FT(e^{j2\pi f_c t}x(t))\right)\\
  &=e^{j2\pi f_c t}x(t)
  \end{aligned}
  \end{equation}
  $$
  这里其实只考虑了单边谱，因此计算增益的时候应该 $\times 2$

2. 低通滤波器的设计

    使用冲击响应不变法，对 $r(t)$ 进行采样，采样间隔为T，即：
    $$
    \begin{equation}
    \begin{aligned}
    r_d(i) = r(iT)= (iT)^{n-1}e^{-2\pi bTi}=T^{n-1}i^{n-1}e^{-2\pi bTi}\label{eq1}
    \end{aligned}
    \end{equation}
    $$
    令 $k=e^{-2\pi bT}$ ，上式就可以简化为
    $$
    \begin{equation}
    \begin{aligned}
    r_d(i)=T^{n-1}i^{n-1}k^{i}
    \end{aligned}
    \end{equation}
    $$

    因为
    $$
    \begin{equation}
    \begin{aligned}
    Z(if(i))&=\sum{if(i)z^{-i}}=-z\frac{\partial\sum{f(i)z^{-i}}}{\partial z}=z^{-1}\frac{\partial\sum{f(i)(z^{-1})^i}}{\partial z^{-1}}\\
    \end{aligned}
    \end{equation}
    $$

    $$
    \begin{equation}
    \begin{aligned}
    Z(k^{i})&=\frac{1}{1-kz^{-1}}
    \end{aligned}
    \end{equation}
    $$
    所以
    $$
    \begin{equation}
    \begin{aligned}
    Z(ik^i)&=z^{-1}\frac{\partial Z(k^{i})}{\partial z^{-1}}=z^{-1}\frac{\partial \frac{1}{1-kz^{-1}}}{\partial z^{-1}}=\frac{kz^{-1}}{(1-kz^{-1})^2}\\

    Z(i^2k^i)&=z^{-1}\frac{\partial Z(iz^{i})}{\partial z^{-1}}=z^{-1}\frac{\frac{kz^{-1}}{(1-kz^{-1})^2}}{\partial z^{-1}}=z^{-1}\left[\frac{k}{(1-kz^{-1})^2}+\frac{2k^2z^{-1}}{(1-kz^{-1})^3}\right]\\

    Z(i^3k^i)&=z^{-1}\frac{\partial Z(i^2k^{i})}{\partial z^{-1}}=z^{-1}\frac{z^{-1}\left[\frac{k}{(1-kz^{-1})^2}+\frac{2k^2z^{-1}}{(1-kz^{-1})^3}\right]}{\partial z^{-1}}\\
    &=z^{-1}\left[\frac{k}{(1-kz^{-1})^2}+\frac{2k^2z^{-1}}{(1-kz^{-1})^3}\right]+z^{-2}\left[\frac{2k^2}{(1-kz^{-1})^3}+\frac{2k^2}{(1-kz)^3}+\frac{6k^3z^{-1}}{(1-kz^{-1})^4}\right]\\
    &=\frac{z^{-1}k(1-kz^{-1})^2+2k^2z^{-2}(1-kz^{-1})+4k^2z^{-2}(1-kz^{-1})+6k^3z^3}{(1-kz^{-1})^4}\\
    &=\frac{kz^{-1}+k^3z^{-3}-2k^2z^{-2}+2k^2z^{-2}-2k^3z^{-3}+4k^2z^{-2}-4k^3z^{-3}+6k^3z^{-3}}{(1-kz^{-1})^4}\\
    &=\frac{k^3z^{-3}+4k^2z^{-2}+kz^{-1}}{(1-kz^{-1})^4}\\
    &=\frac{kz^{-1}(1+4kz^{-1}+k^2z^{-2})}{(1-2kz^{-1}+k^2z^{-2})^2}\\
    &=\frac{kz^{-1}(1+4kz^{-1}+k^2z^{-2})}{1-4kz^{-1}+2k^2z^{-2}+4k^2z^{-2}-4k^3z^{-3}+k^4z^{-4}}\\
    &=\frac{kz^{-1}(1+4kz^{-1}+k^2z^{-2})}{1-4kz^{-1}+6k^2z^{-2}-4k^3z^{-3}+k^4z^{-4}}\\
    \end{aligned}
    \end{equation}
    $$
    因此
    $$
    \begin{equation}
    \begin{aligned}
    Z(r_d(i)) &= Z\left(\frac{(n-1)!}{(2\pi b)^n}T^3*i^3k^i\right)=\frac{6T^3}{(2\pi b)^4}Z(i^3k^i)\\
    &=\frac{6T^3k}{(2\pi b)^4}\frac{z^{-1}(1+4kz^{-1}+k^2z^{-2})}{1-4kz^{-1}+6k^2z^{-2}-4k^3z^{-3}+k^4z^{-4}}\\
    \end{aligned}
    \end{equation}
    $$
    令 $c=\frac{6T^3k}{(2\pi b)^4}$ ,则上式子可重写为
    $$
    \begin{equation}
    \begin{aligned}
    Z(r_d(i))&=c\frac{z^{-1}(1+4kz^{-1}+k^2z^{-2})}{1-4kz^{-1}+6k^2z^{-2}-4k^3z^{-3}+k^4z^{-4}}\\
    \end{aligned}
    \end{equation}
    $$
    最终得到Gammatone滤波器的公式
    $$
    \begin{equation}
    y(n)=c\left[\underbrace{x(n-1)+4kx(n-2)+k^2x(n-3)}_\text{x part}+\underbrace{4ky(n-1)-6k^2y(n-2)+4k^3y(n-3)-k^4y(n-4)}_\text{y part}\right]
    \begin{aligned}
    \end{aligned}
    \end{equation}
    $$

  3. 滤波结果升频
  $$
  \begin{equation}
  \begin{aligned}
  y'(t)=Real\left(IFT(Y(f+f_c))\right)&=Real\left(IFT\left(\int{y(t)e^{-j2\pi(f+f_c)t}dt}\right)\right)\\
  &=Real\left(IFT\left(FT(e^{-j2\pi f_c t}y(t))\right)\right)\\
  &=Real\left(e^{-j2\pi f_c t}y(t)\right)
  \end{aligned}
  \end{equation}
  $$
  升频后信号的实部即为最终结果。

  ## 滤波器中心频率增益归一
  因为实现过程中将带通滤波器转换为低通滤波器，因此只需要对低通滤波器0频处的增益归一为1/2即可。
  因此，归一化系数应该是
  $$
  \begin{equation}
  \begin{aligned}
  scale &= \frac{1}{Z(r_d(i)|_{z=1})/2}=\frac{1}{c\frac{z^{-1}(1+4kz^{-1}+k^2z^{-2})}{1-4kz^{-1}+6k^2z^{-2}-4k^3z^{-3}+k^4z^{-4}}|_{z=1}/2}\\
  &= \frac{(1-k)^4}{c*(1+4k+k^2)}
  \end{aligned}
  \end{equation}
  $$

  ### 误差
  低通滤波器经移频之后，左右两部分可能存在overlap，从而使得带通滤波器中心频率处的增益略大于低通滤波器0频处的增益。误差系数为
  $$
  \begin{equation}
  \begin{aligned}
  \frac{\sqrt{(16Q^4-24Q^2+2)^2+(8Q-32Q^3)^2}}{\sqrt{(16Q^4-24Q^2+1)^2+(8Q-32Q^3)^2}} \approx 1
  \end{aligned}
  \end{equation}
  $$
