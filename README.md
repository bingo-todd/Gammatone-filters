# Gammatone-filters
Python implementation of all-pole Gammatone filter.
The filtering part of code is written in C.
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



[^Holdsworth1988]: Holdsworth, John, Roy Patterson, and Ian Nimmo-Smith. Implementing a GammaTone Filter Bank

[^Brown1994]: G. J. Brown and M. P. Cooke (1994) Computational auditory scene analysis. Computer Speech and Language, 8, pp. 297-336
