# Gammatone-filters
Python implementation of all-pole Gammatone filter.
The filtering part of code is written in C.

Gammatone滤波器冲击响应（Impulse response,IR）：

$$
g(t) = \frac{at^{n-1}\cos(2\pi f_ct+\phi)}{e^{2\pi b t}}
$$

其中:

- $f_c$：中心频率
- $b$ ：带宽，$1.019*ERB(f_c)$


## 中心频率处的增益和相移

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
所以有
$$
\begin{equation}
\begin{aligned}
G(f)&=a \times R(f)*S(f)\\
&=a \times e^{j\phi}\frac{(n-1)!}{(2\pi b)^n}\frac{1}{(1+j(f-f_c)/b)^n}+ae^{-j\phi}\frac{(n-1)!}{(2\pi b)^n}\frac{1}{(1+j(f+f_c)/b)^n}\\
&=a\frac{(n-1)!}{(2\pi b)^n}\left[e^{j\phi}\left(\frac{1}{(1+j(f-f_c)/b)}\right)^n+e^{-j\phi}\left(\frac{1}{(1+j(f+f_c)/b)}\right)^n\right]
\end{aligned}
\end{equation}
$$
对于中心频率处，有
$$
\begin{equation}
\begin{aligned}
\left.G(f)\right|_{f=f_c}&=\left.a\frac{(n-1)!}{(2\pi b)^n}\left[e^{j\phi}\left(\frac{1}{(1+j(f-f_c)/b)}\right)^n+e^{-j\phi}\left(\frac{1}{(1+j(f+f_c)/b)}\right)^n\right]\right|_{f=f_c}\\
&=a\frac{(n-1)!}{(2\pi b)^n}\left[e^{j\phi}+e^{-j\phi}\frac{1}{(1+2jf_c/b)^n}\right]\\
&=a\frac{(n-1)!}{(2\pi b)^n}\left[e^{j\phi}+e^{-j\phi}\frac{1}{(1+2jQ)^n}\right]
\end{aligned}
\end{equation}
$$
通常 $cos(2\pi f_c t+\phi)$ 中的起始相位 $\phi$ 为0，即：
$$
\begin{equation}
\begin{aligned}
Gain(f=f_c)&=\frac{(n-1)!}{(2\pi b)^n }\left[1+\frac{1}{(1+j2Q)^n}\right]\\
&=\frac{6}{(2\pi b)^4}\left[1+\frac{1}{(1+j2f/b)^4}\right]\\
&=\frac{6}{(2\pi b)^4}\left[1+\frac{1}{(1-4Q^2+4jQ)^2}\right]\\
&=\frac{6}{(2\pi b)^4}\left[1+\frac{1}{1-8Q^2+16Q^4-16Q^2+2(1-4Q^2)4jQ}\right]\\
&=\frac{6}{(2\pi b)^4}\left[1+\frac{1}{16Q^4-24Q^2+1+8jQ(1-4Q^2)}\right]\\
&=\frac{6}{(2\pi b)^4}\left[\frac{16Q^4-24Q^2+2+8jQ(1-4Q^2)}{16Q^4-24Q^2+1+8jQ(1-4Q^2)}\right]\\
&=\frac{3}{(2\pi b)^4}\frac{r_1e^{\phi_1}}{r_2e^{\phi_2}}\\
\end{aligned}
\end{equation}
$$
其中
$$
\begin{equation}
\begin{aligned}
\begin{cases}
r_1 = \sqrt{(16Q^4-24Q^2+2)^2+(8Q-32Q^3)^2}\\
\phi_1 = \arctan{\frac{8Q-32Q^3}{16Q^4-24Q^2+2}}\\
r_2 = \sqrt{(16Q^4-24Q^2+1)^2+(8Q-32Q^3)^2}\\
\phi_2 = \arctan{\frac{8Q-32Q^3}{16Q^4-24Q^2+1}}
\end{cases}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{aligned}
Gain_{f_c} = \frac{3}{(2\pi b)^4}\frac{\sqrt{(16Q^4-24Q^2+2)^2+(8Q-32Q^3)^2}}{\sqrt{(16Q^4-24Q^2+1)^2+(8Q-32Q^3)^2}}\\
\phi_{f_c} = \arctan{\frac{8Q-32Q^3}{16Q^4-24Q^2+2}}-\arctan{\frac{8Q-32Q^3}{16Q^4-24Q^2+1}}\\
\end{aligned}
\end{equation}
$$


## Basic ideas [^Holdsworth1988]
Gammatone filter can be regarded as low-pass filter with frequency shitfted by fc. Now, equalently, we can first shift input signal by -fc and filter it with a lowpass filter, finally shift the frequency by fc.
![diagram](images/diagram.png)!

Details, see [README.pdf](README.pdf), currently written in Chinese, but most part are math equations.

## Example

  ```Python
  gtf = APGTF(fs=44100,low_cf=80,high_cf=5000,N_band=32)
  x_filtered = gtf.filter_c(x,is_aligned=0)# not aligned
  ```

### Delays and gains at cfs
  Take filter with cf=4kHz for example, the gain and phase spectrum
  ![delays_gains](images/filter_spectrum.png)

  For delays and gains at center frequencies
  ![delay_gain.png](images/delay_gain.png)
  Basically, the phase delay at center frequency approximates 0.

## Impulse response of Gammatone filters
- Max-amplitude normalized on all bands

  ![ir.png](images/ir.png)

- Gain normalization

  ![ir_norm](images/ir_norm.png)

- Phase compensation

  Phase compensation is actually to align the peaks of all filter impulse response[^Brown1994].

  ![ir_norm_aligned](images/ir_norm_aligned.png)

  <!-- Next, I want to make summary about signal recovery after filtered by Gammatone filters.[Flag] -->

### About efficiency

code
```shell
 $ python GTF.py efficiency

 time consumed(s)
    c         :0.39
    python    :36.23
```


<table>
<tr> <td><img src='test.png'></td> </tr>
<tr> <td><img src='test_phi0.png'></td> </tr> 
</table>


[^Holdsworth1988]: Holdsworth, John, Roy Patterson, and Ian Nimmo-Smith. Implementing a GammaTone Filter Bank

[^Brown1994]: G. J. Brown and M. P. Cooke (1994) Computational auditory scene analysis. Computer Speech and Language, 8, pp. 297-336
