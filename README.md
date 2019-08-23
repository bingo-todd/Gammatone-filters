# Gammatone-filters
Python implementation of all-pole Gammatone filter. 
The filtering part of code is written in C.

## Basic ideas:

Gammatone filter can be regarded as lowpass filter with frequency shitfted by fc. Now, equalently, we can first shift input signal by -fc and filter it with a lowpass filter, finally shift the frequency by fc.

Algorithm details, see [README.pdf](README.pdf), currently written in Chinese,but most part are math equations.

## Impulse response of Gammatone filters
- No phase compensation
![ir_not_aligned](example/ir_not_aligned.png)
- Phase compensation 
![ir_aligned](example/ir_aligned.png)

Next, I want to make summary about signal recovery after filtered by Gammatone filters. ![flag](example/flag.png)
