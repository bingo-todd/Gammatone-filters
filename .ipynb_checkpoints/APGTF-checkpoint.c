//
//  APGTF.c
//  all-pole gammatone filters
//
//  Created by songtao on 2019/8/17.
//  Copyright © 2019 songtao. All rights reserved.
//

#include <stdio.h>
#include<math.h>
#include <stdlib.h>
#include <string.h>
#define MIN_VALUE 1e-20
#define pi 3.14159265358979323846


/*

Args:
    double*x:
    int x_len
    int fs
    double*cfs
    double*bs
    int N_band
    int is_aligned
    int aligned_delay: 
*/
double* APGTF(double*x,int x_len,
              int fs, double*cfs, double*bs, int N_band, 
              int is_aligned, int aligned_delay){
    
    double* y = (double*)malloc(sizeof(double)*x_len*N_band);
    double* gains = (double*)malloc(sizeof(double)*N_band);
    int* delay_array = (int*)malloc(sizeof(int)*N_band);
    int max_delay=0,delay_band=0;
    
    // vars to shift frequency
    double freq_shiftor[2],freq_shiftor_pre[2];
    double cos_step, sin_step;
    
    // vars related to filters
    double tpt = 2*pi*(1.0/(double)fs); // constant
    double p[5][2]; // IIR filter status
    double a[5]; // IIR filter coeffs
    double u0[2]; // FIR filter output
    double b[4]; // FIR filter coeffs
    double k; // constants for each band
    
    int band_i, sample_i, order;
    
    if(is_aligned == 1){
        max_delay = 0;
        for(band_i=0;band_i<N_band;band_i++){
            delay_array[band_i] = round(3.0/(2.0*pi*bs[band_i])*fs);
            if(delay_array[band_i]>max_delay){
                max_delay = delay_array[band_i];
            }
        }
        if(aligned_delay<0){
            aligned_delay = max_delay;
        }
    }
    else{
        for(band_i=0;band_i<N_band;band_i++){
            delay_array[band_i]=0;
        }
        aligned_delay = 0;
    }
    
    memset(y,0,sizeof(double)*x_len*N_band);
    for(band_i=0;band_i<N_band;band_i++){
        // initiation of filter states
        memset(p,0,sizeof(double)*10);
    
        k = exp(-tpt*bs[band_i]);
        // filter coefficients
        a[0] = 1; a[1] = 4*k; a[2] = -6*k*k; a[3]=4*pow(k,3); a[4]=-pow(k, 4);
        b[0] = 1; b[1] = 1;   b[2] = 4*k;    b[3] = k*k;
        
        // filter amp gain normalized
        gains[band_i] = pow(1-k,4)/(1+4*k+k*k)*2;
        
        // computation acceleration
        // convert cos(\phi1+\ph2) and sin(\phi1+\phi2) into mutliplication and summation 
        cos_step = cos(tpt*cfs[band_i]); sin_step = sin(tpt*cfs[band_i]);
        freq_shiftor_pre[0] = cos_step;
        freq_shiftor_pre[1] = -sin_step;
        delay_band = delay_array[band_i]-aligned_delay;
        for(sample_i=0; sample_i<x_len+delay_band; sample_i++){
            freq_shiftor[0] = freq_shiftor_pre[0]*cos_step - freq_shiftor_pre[1]*sin_step;
            freq_shiftor[1] = freq_shiftor_pre[1]*cos_step + freq_shiftor_pre[0]*sin_step;
            freq_shiftor_pre[0] = freq_shiftor[0];
            freq_shiftor_pre[1] = freq_shiftor[1];
            
            // denominator part of filter equation
            // equivalent to add zeros to the end of x
            if(sample_i>=x_len){
                p[0][0] = 0;
                p[0][1] = 0;
            }
            else{
                p[0][0] = x[sample_i]*freq_shiftor[0];
                p[0][1] = x[sample_i]*freq_shiftor[1];
            }
            
            for(order=1;order<=4;order++){
                p[0][0] = p[0][0]+a[order]*p[order][0];
                p[0][1] = p[0][1]+a[order]*p[order][1];
            }
            
            // numerator part of filter equation
            u0[0] = 0; u0[1] = 0;
            for(order=1;order<=3;order++){
                u0[0] = u0[0]+b[order]*p[order][0];
                u0[1] = u0[1]+b[order]*p[order][1];
            }

            // final output = real part of filte result
            if(sample_i>=delay_band){
                y[band_i*x_len+sample_i-delay_band] = gains[band_i]*(u0[0]*freq_shiftor[0] +
                                                                     u0[1]*freq_shiftor[1]);
            }
            // update filter states
            for(order=4;order>=1;order--){
                p[order][0] = p[order-1][0];
                p[order][1] = p[order-1][1];
            }
        }
    }
    return y;
}
