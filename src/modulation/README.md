# Modulation

This directory contains waveform generators and modulation utilities used by RF-EMT. The table below summarizes each source and header file currently in this directory. The `_f64` and `_fp16` variants provide implementations for the indicated floating-point precision; matching headers declare their interfaces.

| File | Summary |
| --- | --- |
| `GMS_IQ_random_bitstream.cpp` | Implements generation of random I/Q bitstream data. |
| `GMS_IQ_random_bitstream.h` | Declares random I/Q bitstream generation interfaces. |
| `GMS_IQ_random_bitstream_f64.cpp` | Implements random I/Q bitstream generation using double-precision data. |
| `GMS_IQ_random_bitstream_f64.h` | Declares the double-precision random I/Q bitstream interfaces. |
| `GMS_am_bb_cmplx_cos_sequence.cpp` | Implements complex baseband amplitude-modulated cosine sequences. |
| `GMS_am_bb_cmplx_cos_sequence.h` | Declares complex baseband AM cosine-sequence interfaces. |
| `GMS_am_bb_cmplx_cos_signal.cpp` | Implements complex baseband AM cosine signals. |
| `GMS_am_bb_cmplx_cos_signal.h` | Declares complex baseband AM cosine-signal interfaces. |
| `GMS_am_bb_cmplx_cos_signal_f64.cpp` | Implements double-precision complex baseband AM cosine signals. |
| `GMS_am_bb_cmplx_cos_signal_f64.h` | Declares double-precision complex baseband AM cosine-signal interfaces. |
| `GMS_am_bb_cmplx_cos_signal_fp16.cpp` | Implements half-precision complex baseband AM cosine signals. |
| `GMS_am_bb_cmplx_cos_signal_fp16.h` | Declares half-precision complex baseband AM cosine-signal interfaces. |
| `GMS_am_bb_cmplx_cosinc_sequence.cpp` | Implements complex baseband AM cosine-increment sequences. |
| `GMS_am_bb_cmplx_cosinc_sequence.h` | Declares complex baseband AM cosine-increment sequence interfaces. |
| `GMS_am_bb_cmplx_cosinc_signal.cpp` | Implements complex baseband AM cosine-increment signals. |
| `GMS_am_bb_cmplx_cosinc_signal.h` | Declares complex baseband AM cosine-increment signal interfaces. |
| `GMS_am_bb_cmplx_cosinc_signal_f64.cpp` | Implements double-precision complex baseband AM cosine-increment signals. |
| `GMS_am_bb_cmplx_cosinc_signal_f64.h` | Declares double-precision complex baseband AM cosine-increment signal interfaces. |
| `GMS_am_bb_cmplx_sinc_sequence.cpp` | Implements complex baseband AM sinc sequences. |
| `GMS_am_bb_cmplx_sinc_sequence.h` | Declares complex baseband AM sinc-sequence interfaces. |
| `GMS_am_bb_cmplx_sinc_signal.cpp` | Implements complex baseband AM sinc signals. |
| `GMS_am_bb_cmplx_sinc_signal.h` | Declares complex baseband AM sinc-signal interfaces. |
| `GMS_am_bb_cmplx_sinc_signal_f64.cpp` | Implements double-precision complex baseband AM sinc signals. |
| `GMS_am_bb_cmplx_sinc_signal_f64.h` | Declares double-precision complex baseband AM sinc-signal interfaces. |
| `GMS_am_bb_cmplx_sine_signal.cpp` | Implements complex baseband AM sine signals. |
| `GMS_am_bb_cmplx_sine_signal.h` | Declares complex baseband AM sine-signal interfaces. |
| `GMS_am_bb_cmplx_sine_signal_f64.cpp` | Implements double-precision complex baseband AM sine signals. |
| `GMS_am_bb_cmplx_sine_signal_f64.h` | Declares double-precision complex baseband AM sine-signal interfaces. |
| `GMS_am_bb_cmplx_sine_signal_fp16.cpp` | Implements half-precision complex baseband AM sine signals. |
| `GMS_am_bb_cmplx_sine_signal_fp16.h` | Declares half-precision complex baseband AM sine-signal interfaces. |
| `GMS_am_bb_cmplx_trapez_signal.cpp` | Implements complex baseband AM trapezoidal signals. |
| `GMS_am_bb_cmplx_trapez_signal.h` | Declares complex baseband AM trapezoidal-signal interfaces. |
| `GMS_am_bb_cmplx_trapez_signal_f64.cpp` | Implements double-precision complex baseband AM trapezoidal signals. |
| `GMS_am_bb_cmplx_trapez_signal_f64.h` | Declares double-precision complex baseband AM trapezoidal-signal interfaces. |
| `GMS_am_bb_cmplx_trapez_signal_fp16.cpp` | Implements half-precision complex baseband AM trapezoidal signals. |
| `GMS_am_bb_cmplx_trapez_signal_fp16.h` | Declares half-precision complex baseband AM trapezoidal-signal interfaces. |
| `GMS_am_bb_cosine_sequence.cpp` | Implements baseband AM cosine sequences. |
| `GMS_am_bb_cosine_sequence.h` | Declares baseband AM cosine-sequence interfaces. |
| `GMS_am_bb_cosine_signal.cpp` | Implements baseband AM cosine signals. |
| `GMS_am_bb_cosine_signal.h` | Declares baseband AM cosine-signal interfaces. |
| `GMS_am_bb_cosine_signal_f64.cpp` | Implements double-precision baseband AM cosine signals. |
| `GMS_am_bb_cosine_signal_f64.h` | Declares double-precision baseband AM cosine-signal interfaces. |
| `GMS_am_bb_cosine_signal_fp16.cpp` | Implements half-precision baseband AM cosine signals. |
| `GMS_am_bb_cosine_signal_fp16.h` | Declares half-precision baseband AM cosine-signal interfaces. |
| `GMS_am_bb_sine_sequence.cpp` | Implements baseband AM sine sequences. |
| `GMS_am_bb_sine_sequence.h` | Declares baseband AM sine-sequence interfaces. |
| `GMS_am_bb_sine_signal.cpp` | Implements baseband AM sine signals. |
| `GMS_am_bb_sine_signal.h` | Declares baseband AM sine-signal interfaces. |
| `GMS_am_bb_sine_signal_f64.cpp` | Implements double-precision baseband AM sine signals. |
| `GMS_am_bb_sine_signal_f64.h` | Declares double-precision baseband AM sine-signal interfaces. |
| `GMS_am_bb_sine_signal_fp16.cpp` | Implements half-precision baseband AM sine signals. |
| `GMS_am_bb_sine_signal_fp16.h` | Declares half-precision baseband AM sine-signal interfaces. |
| `GMS_am_bb_square_signal.cpp` | Implements baseband AM square signals. |
| `GMS_am_bb_square_signal.h` | Declares baseband AM square-signal interfaces. |
| `GMS_am_bb_square_signal_f64.cpp` | Implements double-precision baseband AM square signals. |
| `GMS_am_bb_square_signal_f64.h` | Declares double-precision baseband AM square-signal interfaces. |
| `GMS_am_bb_square_signal_fp16.cpp` | Implements half-precision baseband AM square signals. |
| `GMS_am_bb_square_signal_fp16.h` | Declares half-precision baseband AM square-signal interfaces. |
| `GMS_am_bb_trapez_signal.cpp` | Implements baseband AM trapezoidal signals. |
| `GMS_am_bb_trapez_signal.h` | Declares baseband AM trapezoidal-signal interfaces. |
| `GMS_am_bb_trapez_signal_f64.cpp` | Implements double-precision baseband AM trapezoidal signals. |
| `GMS_am_bb_trapez_signal_f64.h` | Declares double-precision baseband AM trapezoidal-signal interfaces. |
| `GMS_am_bb_trapez_signal_fp16.cpp` | Implements half-precision baseband AM trapezoidal signals. |
| `GMS_am_bb_trapez_signal_fp16.h` | Declares half-precision baseband AM trapezoidal-signal interfaces. |
| `GMS_amplitude_modulator.cpp` | Implements amplitude-modulation operations. |
| `GMS_amplitude_modulator.h` | Declares amplitude-modulator interfaces. |
| `GMS_bit_reverse_simd.h` | Provides SIMD helpers for reversing bit order. |
| `GMS_cmplx_trapezw_env.cpp` | Implements a complex trapezoidal waveform envelope. |
| `GMS_cmplx_trapezw_env.h` | Declares complex trapezoidal-envelope interfaces. |
| `GMS_cpfsk_modulator.h` | Declares continuous-phase frequency-shift keying (CPFSK) modulation interfaces. |
| `GMS_cpfsk_signal.cpp` | Implements CPFSK signal generation. |
| `GMS_cpfsk_signal.h` | Declares CPFSK signal-generation interfaces. |
| `GMS_cpm_pulse_shapers.cpp` | Implements pulse-shaping functions for continuous-phase modulation (CPM). |
| `GMS_cpm_pulse_shapers.h` | Declares CPM pulse-shaping interfaces. |
| `GMS_dsp_simd_fix_float_kernels.h` | Provides SIMD DSP kernels for fixed-point and floating-point data. |
| `GMS_iq_rectw_bitstream_vsequence.cpp` | Implements I/Q rectangular-wave bitstream vector sequences. |
| `GMS_iq_rectw_bitstream_vsequence.h` | Declares I/Q rectangular-wave bitstream vector-sequence interfaces. |
| `GMS_iq_rectwave_bitstream.cpp` | Implements I/Q rectangular-wave bitstream generation. |
| `GMS_iq_rectwave_bitstream.h` | Declares I/Q rectangular-wave bitstream interfaces. |
| `GMS_ofdm_signal_modulator.h` | Declares an orthogonal frequency-division multiplexing (OFDM) signal modulator. |
| `GMS_pm_bb_cmplx_cos_signal.cpp` | Implements complex baseband phase-modulated cosine signals. |
| `GMS_pm_bb_cmplx_cos_signal.h` | Declares complex baseband PM cosine-signal interfaces. |
| `GMS_pm_bb_cmplx_cos_signal_f64.cpp` | Implements double-precision complex baseband PM cosine signals. |
| `GMS_pm_bb_cmplx_cos_signal_f64.h` | Declares double-precision complex baseband PM cosine-signal interfaces. |
| `GMS_pm_bb_cmplx_cosinc_signal.cpp` | Implements complex baseband phase-modulated cosine-increment signals. |
| `GMS_pm_bb_cmplx_cosinc_signal.h` | Declares complex baseband PM cosine-increment signal interfaces. |
| `GMS_pm_bb_cmplx_sinc_signal.cpp` | Implements complex baseband phase-modulated sinc signals. |
| `GMS_pm_bb_cmplx_sinc_signal.h` | Declares complex baseband PM sinc-signal interfaces. |
| `GMS_pm_bb_cmplx_sine_signal.cpp` | Implements complex baseband phase-modulated sine signals. |
| `GMS_pm_bb_cmplx_sine_signal.h` | Declares complex baseband PM sine-signal interfaces. |
| `GMS_pm_bb_cmplx_trapez_signal.cpp` | Implements complex baseband phase-modulated trapezoidal signals. |
| `GMS_pm_bb_cmplx_trapez_signal.h` | Declares complex baseband PM trapezoidal-signal interfaces. |
| `GMS_pm_bb_cmplx_trapez_signal_fp16.cpp` | Implements half-precision complex baseband PM trapezoidal signals. |
| `GMS_pm_bb_cmplx_trapez_signal_fp16.h` | Declares half-precision complex baseband PM trapezoidal-signal interfaces. |
| `GMS_rand_pm_gaussian_signal_f32.cpp` | Implements random phase-modulated Gaussian signals using single-precision data. |
| `GMS_rand_pm_gaussian_signal_f32.h` | Declares single-precision random PM Gaussian-signal interfaces. |
| `GMS_rand_pm_gaussian_signal_f64.cpp` | Implements random phase-modulated Gaussian signals using double-precision data. |
| `GMS_rand_pm_gaussian_signal_f64.h` | Declares double-precision random PM Gaussian-signal interfaces. |
| `GMS_rectangular_waveform.cpp` | Implements rectangular waveform generation. |
| `GMS_rectangular_waveform.h` | Declares rectangular waveform interfaces. |
| `GMS_rectangular_waveform_fp16.cpp` | Implements rectangular waveform generation using half-precision data. |
| `GMS_rectangular_waveform_fp16.h` | Declares half-precision rectangular waveform interfaces. |
| `GMS_sawtooth_waveform.cpp` | Implements sawtooth waveform generation. |
| `GMS_sawtooth_waveform.h` | Declares sawtooth waveform interfaces. |
| `GMS_sinusoidal_fsk.cpp` | Implements sinusoidal frequency-shift keying (FSK) modulation. |
| `GMS_sinusoidal_fsk.h` | Declares sinusoidal FSK modulation interfaces. |
| `GMS_sinusoidal_weighted_oqpsk.cpp` | Implements sinusoidally weighted offset QPSK modulation. |
| `GMS_sinusoidal_weighted_oqpsk.h` | Declares sinusoidally weighted offset QPSK interfaces. |
| `GMS_square_waveform.cpp` | Implements square waveform generation. |
| `GMS_square_waveform.h` | Declares square waveform interfaces. |
| `GMS_trapezoid_waveform.cpp` | Implements trapezoidal waveform generation. |
| `GMS_trapezoid_waveform.h` | Declares trapezoidal waveform interfaces. |
| `GMS_trapezoid_waveform_fp16.cpp` | Implements trapezoidal waveform generation using half-precision data. |
| `GMS_trapezoid_waveform_fp16.h` | Declares half-precision trapezoidal waveform interfaces. |
| `GMS_triangle_waveform.cpp` | Implements triangle waveform generation. |
| `GMS_triangle_waveform.h` | Declares triangle waveform interfaces. |
| `GMS_white_gauss_noise.cpp` | Implements white Gaussian noise generation. |
| `GMS_white_gauss_noise.h` | Declares white Gaussian noise generation interfaces. |
| `GMS_white_gauss_noise.inl` | Provides inline implementation details for white Gaussian noise generation. |
