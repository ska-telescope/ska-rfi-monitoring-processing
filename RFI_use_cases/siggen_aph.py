"""
    Work contributed to the SKA Telescope project (http://www.skatelescope.org) by:
    SKA South Africa (http://www.ska.ac.za)
    
    A library for generating signals for simulations.
    Note: the reference impedance in this module that relates amplitude to power, is Z0=1 by default.
    
    @author: aph@ska.ac.za
    @last_modified: 16-03-2018
    
"""

from __future__ import division
import copy as cp
import numpy as np
import scipy.signal

try:
    from matplotlib import mlab, pyplot
    import atexit # Before exiting make sure to show plots and pause
    def _show_at_exit_():
        fignums = pyplot.get_fignums()
        if (len(fignums) > 0):
            pyplot.show()
    atexit.register(_show_at_exit_)
except:
    print("Failed to import matplotlib, won't be able to show()!")


kB = 1.38e-23 # Boltzmann's constant
MHz = 1e6 # Convenient constants
Z0 = 1 # Reference impedance used to relate amplitude to power


def signal_power(signal):
    """
        Uses globally defined Z0 to relate amplitude and power.
        @param signal: time domain input signal series.
        @return: the total integrated power of the signal [W] (linear scale)
    """
    global Z0
    return np.std(signal)**2/Z0

def scale_signal(signal, P_signal):
    """
        Scales the amplitude of a signal for the specified power (w.r.t Z0).
        @param signal: time domain signal series.
        @param P_signal: the total integrated power in the particular component [W] (linear scale).
        @return: scaled copy of signal.
    """
    P_signal0 = signal_power(signal)
    return signal*(P_signal/P_signal0)**.5


def band_limit(signal, samplerate, f_pass, f_stop, pass_dB, stop_dB, ftype="cheby1", TFT='sos'):
    """
        Band pass filtering using 'scipy.signal.sosfiltfilt()' or 'scipy.signal.filtfilt()' to yield a net zero
        phase response. The signal is zero padded by len(signal)-1 with "odd" symmetry, to minimize the filter's
        transient effects.
        Note: TFT='sos' -> cascaded second order sections are numerically more robust for high order filters than 'b/a'.
        @param samplerate: sample frequency of "s_sampled"
        @param f_pass, f_stop: pass band & stop band frequencies in the analogue domain (same units as "samplerate")
        @param pass_dB: filter permissible ripple (pass-band or stop-band or both) [dB].
        @param stop_dB: filter minimum attenuation at stopband frequencies [dB]
        @param ftype: 'cheby1' has fast stopband transition with monotonic response, but passband ripple.
                      'cheby2' has monotonic/flat passband with slower stopband transition with stopband ripple.
                      'ellip' has fastest stopband transition with ripple in both pass and stopbands and usually
                       has lowest order of all types for similar characteristics.
                      (default 'cheby1').
        @param TFT: 'sos' for cascaded 2nd order realisation, or 'ba' for generic transfer function (default 'sos' since
                    'ba' sometimes fails silently if filter order becomes too large!)
        @return: band-pass filtered version of signal.
    """
    assert (max(f_pass) <= samplerate/2.), "Passband frequencies cannot exceed Nyquist rate"
    assert (max(f_stop) <= samplerate/2.), "Stopband frequencies cannot exceed Nyquist rate"
    wp = [2*f_pass[0]/samplerate, 2*f_pass[1]/samplerate]
    ws = [2*f_stop[0]/samplerate, 2*f_stop[1]/samplerate]
    if (wp[0] == 0) and (ws[0] == 0): # It's a low pass filter
        wp = wp[1]
        ws = ws[1]
    elif (wp[1] == 0.5) and (ws[1] == 0.5): # It's a high pass filter
        wp = wp[0]
        ws = ws[0]
    # Attenuation is doubled by filtfilt applying twice
    coeff = scipy.signal.iirdesign(wp, ws, gpass=pass_dB/2., gstop=stop_dB/2., analog=0, ftype=ftype, output=TFT)
    # filtfilt() applies the filter twice, in steady state
    if (TFT=='ba'):
        y = scipy.signal.filtfilt(coeff[0],coeff[1], signal, padlen=len(signal)-1)
    elif (TFT=='sos'):
        y = scipy.signal.sosfiltfilt(coeff, signal, padlen=len(signal)-1)
    return y


def WhiteNoiseSignal(t_sample, dBm_Hz=None, Teq=None):
    """
        @param t_sample: time series for which to generate the signal [sec].
        @param dBm_Hz: power spectral density [dBm/Hz].
        @param Teq: noise equivalent temperature [K], only used if dBm_Hz not given.
        @return: signal series
    """
    samplerate = 1/np.diff(t_sample).mean()
    psd = (10**(dBm_Hz/10.) * 1e-3) if (dBm_Hz is not None) else (kB*Teq) 
    return scale_signal(np.random.randn(len(t_sample)), psd*samplerate/2.)


def PinkNoiseSignal(t_sample, dBm, pass_MHz, stop_MHz, stop_dB):
    """
        @param t_sample: time series for which to generate the signal [sec].
        @param dBm: total integrated power [dBm].
        @param pass_MHz: 0.5dB bandwidth as (f_low,f_high) [MHz].
        @param stop_MHz: stop-band frequencies as (f_low,f_high) [MHz].
        @param stop_dB: attenuation at stop-band frequencies relative to pass-band frequencies [dB].
        @return: signal series.
    """
    samplerate = 1/np.diff(t_sample).mean()
    s = band_limit(WhiteNoiseSignal(t_sample,Teq=1e9), samplerate/MHz, pass_MHz, stop_MHz, pass_dB=3, stop_dB=stop_dB, ftype="butter")
    return scale_signal(s, 10**(dBm/10.)*1e-3)


def PulsedSignal(t_sample, dBm_pk, f_MHz, T_separate, pattern=None):
    """
        A repeating pulse train with first peak at t_sample=T_separate/2 (or overruled by "pattern").
        @param t_sample: time series for which to generate the signal [sec].
        @param dBm_pk: the power of the carrier signal [dBm].
        @param f_MHz: the carrier frequency [MHz].
        @param T_separate: the minimum time interval between pulse crests [sec].
        @param pattern: a binary string where "1" identifies an ON and "0" an OFF pulse, None for random.
        @return: signal series.
    """
    V = np.sin(2*np.pi*f_MHz*MHz*t_sample)
    V = scale_signal(V, 10**(dBm_pk/10.)*1e-3)
    
    if (pattern is not None): # Blank out pulses according to pattern
        pattern = [int(p) for p in pattern]
    else: # Blank out some pulses with random {0,1} symbols
        N_tot = int((t_sample[-1]-t_sample[0])/T_separate)+1 # Total number of symbols (rounded up)
        pattern = np.random.random_integers(0, 1, N_tot) # Rounded no. of symbols to span the time range
        pattern[0] = 1 # Force starting with one
    T_pattern = T_separate*len(pattern)
    pattern = np.tile(pattern, int((np.max(t_sample)-np.min(t_sample))/T_pattern+1)) # Pattern repeats
    samplerate = 1/np.diff(t_sample).mean()
    mask = np.repeat(pattern, int(T_separate*samplerate+0.5))[:len(t_sample)] # Multiple samples per symbol interval
    mask = np.roll(mask, -int(t_sample[0]*samplerate)) # Shift the mask from being aligned with t[0] = 0 to t[0] = t_sample[0]
    
    return V*mask


def GaussianPulsedSignal(t_sample, dBm_pk, f_MHz, T_width, T_separate, pattern=None):
    """
        A repeating pulse train with Gaussian envelope with first peak at t_sample=T_separate/2 (or overruled by "pattern").
        @param t_sample: time series for which to generate the signal [sec].
        @param dBm_pk: the power of the carrier signal [dBm].
        @param f_MHz: the carrier frequency [MHz].
        @param T_width: the pulse half-power width [sec].
        @param T_separate: the minimum time interval between pulse crests [sec].
        @param pattern: a binary string where "1" identifies an ON and "0" an OFF pulse, None for random.
        @return: signal series.
    """
    V = PulsedSignal(t_sample, dBm_pk, f_MHz, T_separate, pattern)
    envelope = np.exp(-2*np.log(4)*((t_sample%T_separate-T_separate/2.)/T_width)**2) # Regular pulse train
    return V*envelope


def symbol_series(t_sample, N_symbols, T_symbol):
    """
        Generates a sampled series of random symbols in non-return-to-zero format.
        @param t_sample: time series for which to generate the signal [sec].
        @param N_symbols: the number of discrete symbols to choose from.
        @param T_symbol: the duration per symbol [sec].
        @return: sampled random symbols, integers from {0..N_symbols-1}
    """
    samplerate = 1/np.diff(t_sample).mean()
    N_tot = int((t_sample[-1]-t_sample[0])/T_symbol)+1 # Total number of symbols (rounded up)
    data = np.random.random_integers(0, N_symbols-1, N_tot) # Rounded no. of symbols to span the time range
    data = np.repeat(data, int(round(T_symbol*samplerate+0.5)))[:len(t_sample)] # Multiple samples per symbol interval
    return data

def ASK_series(t_sample, f_MHz, T_symbol):
    """
        Generates an Amplitude Shift Keying modulated signal based on a random data series.
        @param t_sample: time series for which to generate the signal [sec].
        @param f_MHz: the frequency of the carrier [MHz].
        @param T_symbol: the duration per symbol [sec].
        @return: the sampled ASK signal series.
    """
    # The random series of binary symbols
    data = symbol_series(t_sample, 2, T_symbol)
    # Amplitude-modulate the carrier with the symbols
    samples = data * np.sin(2*np.pi*f_MHz*MHz*t_sample)
    return samples

def PSK_series(t_sample, f_MHz, T_symbol, N_symbols):
    """
        Generates an N-Phase Shift Keying modulated signal based on a random data series.
        @param t_sample: time series for which to generate the signal [sec].
        @param f_MHz: the frequency of the carrier [MHz].
        @param T_symbol: the duration per symbol [sec].
        @param N_symbols: the number of discrete symbols to modulate with.
        @return: the sampled PSK signal series.
    """
    # The random series of symbols
    data = symbol_series(t_sample, N_symbols, T_symbol)
    # Phase-modulate the carrier with the symbols
    samples = np.sin(2*np.pi*f_MHz*MHz*t_sample + 2*np.pi*(data-0.5)/N_symbols)
    return samples

def FSK_series(t_sample, f_MHz, T_symbol, N_symbols, ff=0.1):
    """
        Generates an N-Frequency Shift Keying modulated signal based on a random data series.
        @param t_sample: time series for which to generate the signal [sec].
        @param f_MHz: the frequency of the carrier [MHz].
        @param T_symbol: the duration per symbol [sec].
        @param N_symbols: the number of discrete symbols to modulate with.
        @param ff: the modulation depth as the fractional maximum excursion from "f_MHz" ("BW/2/f_MHz").
        @return: the sampled FSK signal series.
    """
    # The random series of symbols
    data = symbol_series(t_sample, N_symbols, T_symbol)
    # Frequency-modulate the carrier with the symbols
    samples = np.sin(2*np.pi * f_MHz*MHz*(1+ff*(data/(N_symbols-1.)-0.5)) * t_sample)
    return samples

def QAM_series(t_sample, f_MHz, T_symbol, N_symbols):
    """
        Generates an N-Quadrature Amplitude Modulated signal based on a random data series.
        @param t_sample: time series for which to generate the signal [sec].
        @param f_MHz: the frequency of the carrier [MHz].
        @param T_symbol: the duration per symbol [sec].
        @param N_symbols: the number of discrete symbols to modulate with, must be 2N.
        @return: the sampled QAM signal series.
    """
    # Set up the maps of discrete symbols to amplitude & phase
    shell_counts = [(2*s+2)**2 for s in range(int((N_symbols/np.pi)**.5)+1)] # Number of symbols permitted within each shell (up to and including)
    a_map = [np.searchsorted(shell_counts,s+1)+1 for s in range(N_symbols)] # Amplitudes for each symbol, from the shell it's allocated to
    p_map = []
    for c in [shell_counts[0]]+np.diff(shell_counts).tolist(): # Counts per shell
        p_map.extend(np.cumsum([2*np.pi/c]*c)) # Phases for each symbol, computed per shell
    
    # The random series of symbols
    N_tot = int((t_sample[-1]-t_sample[0])/T_symbol)+1 # Total number of symbols (rounded up)
    symbol_data = np.random.random_integers(0, N_symbols-1, N_tot) # Rounded no. of symbols to span the time range
    # Map the symbols to amplitude & phase
    a_data = [a_map[s] for s in symbol_data]
    p_data = [p_map[s] for s in symbol_data]
    
    # Expand to multiple samples per symbol interval
    samplerate = 1/np.diff(t_sample).mean()
    amplitude = np.repeat(a_data, int(round(T_symbol*samplerate+0.5)))[:len(t_sample)] # Amplitude series
    phase = np.repeat(p_data, int(round(T_symbol*samplerate+0.5)))[:len(t_sample)] # Phase series
    
    # Amplitude & phase modulate the carrier with the symbols
    samples = amplitude * np.sin(2*np.pi*f_MHz*MHz*t_sample + phase+np.pi/4.)
    return samples


def _make_TDMA_slot_(t_sample, signal, T_slot, slot_no, N_slots):
    """
        Chops the given signal into the specified TDMA channel.
        @param t_sample: time series for which to generate the signal [sec].
        @param T_slot: time duration of a single slot [sec].
        @param slot_no: the sequence number of this particular time slot, starting at 0.
        @param N_slots: the number of slots in a frame.
        @param winfunc: the time-domain windowing function as "lambda M: winvals" (default Hamming).
        @return: the signal after being limited in time domain.
    """
    # Construct the mask for this slot
    mask = np.zeros_like(signal)
    samplerate = 1/np.diff(t_sample).mean()
    N = int(T_slot*N_slots*samplerate) # Number of time samples in an entire frame
    t_0 = slot_no*T_slot # Time offset for this slot
    n_0 = int(t_0*samplerate) # Offset of this slot in number of time samples
    wlen = N//N_slots
    for F in range(len(t_sample)//N+1): # For each frame
        mask[F*N+n_0:F*N+n_0+wlen] = 1
    # Apply the mask
    signal = signal*mask
    return signal


def package_TDMA_series(t_sample, signals, T_slot, slot_occupancy, f_c=None,BW=None,df_stop=None,stop_dB=None):
    """
        Creates a time series composed of frames which consists of timeslots @ T_slot per timeslot.
        The individual slots of each frame is ON or OFF as per frame_slots. The frame sequences are repeated to fill time.
        Out-of-band suppression is performed after packaging, only if parameters are given.
        @param signals: [signal] time series matching the length of the slot_occupancy pattern.
        @param T_slot: time duration of a single slot [sec].
        @param slot_occupancy: a string consisting of "0" and "1"s, one character for each slot.
        @param f_c, BW, df_stop: frequencies for filtering, in [Hz] (default None).
        @param stop_dB: out-of-band attenuation that's required at f_c+-stop_BW (default None).
        @return: the signal after being limited in time domain.
    """
    V_slots = [_make_TDMA_slot_(t_sample,signals[s],T_slot,s,len(slot_occupancy)) \
                 for s,_sON_ in enumerate(slot_occupancy) if _sON_=="1"]
    V_frame = np.sum(V_slots, axis=0)
    if (f_c and BW and df_stop and stop_dB): # Apply out-of-band attenuation if requested
        samplerate = 1/np.diff(t_sample).mean()
        V_frame = band_limit(V_frame, samplerate, (f_c-BW/2.,f_c+BW/2.), (f_c-df_stop,f_c+df_stop), pass_dB=0.1,stop_dB=stop_dB, ftype='cheby1')
    return V_frame


class sig_generator(object):
    """ A stateless generator to produce signals in blocks. Example of use:
            sig_generator(WhiteNoiseSignal, Teq=13)(t_s)
    """
    def __init__(self, make_signal, **kwargs):
        """
            @param make_signal: a function to generate the signal time series, like lambda time,**kwargs: time_series_amplitude
                    with t_sample: time series for which to generate the signal [sec]. If this is a list of functions, it is assumed
                    that their results must be summed.
            @param kwargs: arguments to be passed to make_signal
        """
        self.make_signal = make_signal
        self.kwargs = kwargs
        self.t_delay = 0
    
    def delay(self, t_delay):
        """ Returns a time-delayed copy of this generator, on top of any existing delay.
            CAUTION: for non-CW signals generated by numpy's random module, the delay is actually an ADVANCE! 
            @param t_delay: the interval to increase the delay of this generator by [sec]
        """
        copy = cp.copy(self)
        try: # If this is a combined generator, push delays down (necessary because of how the random number generator gets delayed)
            copy.make_signal = [copy.make_signal[i].delay(t_delay) for i in range(len(copy.make_signal))]
        except: # It's an atomic generator, accept the delay
            copy.t_delay += t_delay
        return copy
    
    def __call__(self, time):
        """ Generates the (time delayed) signal corresponding to the time series, including delay.
            @return: signal series returned by 'make_signal' function
        """
        # Ensure that all signals based on random sequences are also repeatable
        dt = np.diff(time).mean()
        np.random.seed( (2**31+int(time[0]/dt)) % 2**32) # Add 2**31 to accommodate t[0] being negative
        for z in range(abs(int(self.t_delay/dt))): # ADVANCE the random generator because I can't figure out how to DELAY it.
            np.random.random()
        # Adjust time reference without affecting random repeatability!
        time = time - self.t_delay
        
        # Construct the signal & return
        return np.sum(make_signal(time, **self.kwargs) for make_signal in np.atleast_1d(self.make_signal))
    
    def __add__(self, other_generator):
        """ Returns a new generator object that combines two generator objects by summing their generated signals """
        return sig_generator(make_signal=[self,other_generator])


class pat_generator(sig_generator):
    """ A stateless generator to produce signals in blocks """
    def __init__(self, make_signal, pattern, t_start=None, t_stop=None):
        """
            @param make_signal: a function to generate the signal time series, like lambda time,pattern: time_series_amplitude
                    with t_sample: time series for which to generate the signal [sec].
            @param pattern: either a string representing the sequence of discrete symbols, or an array which can be
                    trivially converted to such a  string.
        """
        sig_generator.__init__(self, make_signal)
        self.master_pattern = "".join(map(str,pattern)) # Iterables become strings 
        self.t_start, self.t_stop = t_start, t_stop
    
    def __call__(self, time):
        """ Generates the (time delayed) signal with pattern matched to the time series.
            @return: signal series returned by 'make_signal' function
        """
        # Generate the sub-interval of pattern to generate signal for
        N_mastersymbols = len(self.master_pattern)
        T_symbol = (self.t_stop-self.t_start)/float(N_mastersymbols)
        N0 = int((time[0]-self.t_start)/T_symbol) % N_mastersymbols # Unwraps to 0..N_mastersymbols
        N_symbols = int((time[-1]-time[0])/T_symbol)
        # This sequence starts at N0, may need to repeat the master pattern if the sequence ends beyond the master pattern
        N_repeat = 1 + int((N0+N_symbols)/float(N_mastersymbols)+0.5)
        master_pattern = self.master_pattern*N_repeat
        pattern = master_pattern[N0:N0+N_symbols+1]
        
        # Construct the signal & return
        self.kwargs["pattern"] = pattern
        return sig_generator.__call__(self, time)


def generate_SSR(dBm_pk, T_repeat=100e-6):
    """
        Random series of Pulse Amplitude Modulated signals, band limited to ~2 MHz bandwidth. First pulse peak at ~0.77 microsec.
        [1] https://en.wikipedia.org/wiki/Secondary_surveillance_radar
        SSR Mode A & C replies (air-to-ground): 12 x 0.45 microsec pulses, separated by 1.45 microsec
        @param dBm_pk: the peak integrated power [dBm].
        @param T_repeat: repeat interval for the pulse train (default 100e-6) [sec].
        @return the signal series.
    """
    T_separate = 1.45e-6
    pattern = "1" + ("10001010011100010100101011100011010110"[np.random.randint(0,15):])
    pattern = pattern[:12] + ("0"*int((T_repeat-18e-6)/T_separate)) # 12 random pulses every T_repeat sec
    def make_signal(t_sample, pattern):
        return GaussianPulsedSignal(t_sample, dBm_pk, f_MHz=1090, T_width=0.45e-6, T_separate=T_separate, pattern=pattern)
    return pat_generator(make_signal, pattern, 0, T_repeat)


def generate_SSR_S(dBm_pk, extended=False, T_repeat=1e-3):
    """
        Random series of SSR-S compliant Pulse Position Modulated signals. First pulse peak at 0.5 microsec.
        [1] https://en.wikipedia.org/wiki/Secondary_surveillance_radar
        [2] http://www.radartutorial.eu/13.ssr/sr24.en.html
        SSR Mode S replies (air-to-ground): 56 x 0.5 microsec pulses, separated by 1.0 microsec
        ADS-B = SSR Mode S "extended" to 112 symbols.
        @param dBm_pk: the peak integrated power [dBm].
        @param extended: False to generate 56 symbol sequence or True to generate 112 (default False).
        @param T_repeat: repeat interval for the pulse train (default 1e-3) [sec].
        @return the signal series.
    """
    T_separate = 1e-6
    N_block = 112 if extended else 56
    encode = lambda ud: ud.replace("u","01").replace("d","10")
    preamble = "dd00dd000000"
    pattern = "udddududduuudddududdududuuuddduududuud"[np.random.randint(0,15):]
    pattern = encode(preamble + (pattern[:14]*int(N_block/14))) + ("00"*(int((T_repeat-8e-6)/T_separate)-N_block))
    def make_signal(t_sample, pattern):
        S = PulsedSignal(t_sample, dBm_pk, f_MHz=1090, T_separate=T_separate/2., pattern=pattern)
        # S: Sidebands suppressed to -10dBC @ +/-1MHz, -30dBC @+/- 7MHz, -40dBC @ +/- 23 MHz, as per regulations, however -50dBC @ +/- 78MHz so add some filtering
        f_sample = 1/np.diff(t_sample).mean()
        return band_limit(S, f_sample/MHz, (1090-23,1090+23), (1090-78,1090+78), pass_dB=3,stop_dB=10, ftype="butter")
    return pat_generator(make_signal, pattern, 0, T_repeat)


def generate_DME(dBm_pk, f_MHz=1115, pattern="10100"):
    """
        Pulse pairs, each 3.5 us wide, separated by 12 us. First pulse peak at 6 us.
        @param dBm_pk: the peak integrated power [dBm].
        @param f_MHz: the centre frequency for the signal [MHz].
        @param pattern: a string consisting of "1" & "0"s, each marking a *pair* of pulses ON & OFF respectively.
        @return the signal series.
    """
    pattern = np.repeat([int(p) for p in pattern], 2) # On & Off symbols occur in pairs
    def make_signal(t_sample, pattern):
        return GaussianPulsedSignal(t_sample, dBm_pk, f_MHz=f_MHz, T_width=3.5e-6, T_separate=12e-6, pattern=pattern)
    return pat_generator(make_signal, pattern, 0, len(pattern)*12e-6)


def generate_RADAR(dBm_pk, f_MHz, T_width=2e-6, T_repeat=2e-3):
    """
        Perfectly repetitive cycle of pulses. First pulse peak at T_width.
        @param dBm_pk: the peak integrated power [dBm].
        @param f_MHz: the centre frequency for the signal [MHz].
        @param T_width, T_repeat: time intervals (default 2microsec & 2millisec) [sec]
        @return the signal series.
    """
    def make_signal(t_sample, pattern):
        return GaussianPulsedSignal(t_sample+(T_repeat/2.-T_width), dBm_pk, f_MHz=f_MHz, T_width=T_width, T_separate=T_repeat, pattern="1")
    return pat_generator(make_signal)


def generate_TDMA(dBm, f_MHz, T_symbol, T_slot, slot_occupancy="1011000", NPSK=4, NQAM=0):
    """
        Narrow band TDMA communications channels at the specified channel frequency.
        The time series is composed of frames which consists of timeslots @ T_slot per timeslot.
        The individual slots of each frame is ON or OFF as per frame_slots. The frame sequences are repeated to fill time.
        
        According to Wikipedia PSK is widely used in existing communication technologies:
        Bluetooth-2 uses 4-PSK (2Mbps) & 8-PSK (3Mbps), Satellite downlinks frequently employ 4-PSK & 8-PSK (HD video, DVB-S2),
        Wireless LAN uses Orthogonal-FDM where each channel is modulated either with 2-PSK (1Mbps) or 4-PSK (2Mbps - 11Mbps).
        Variants of 4-PSK has been adopted for TDMA mobile telephony. Typically for more than 8 symbols per chip, QAM is
        favoured e.g. IMT-2000 (UMTS).
        @param dBm, f_MHz, T_symbol: controls the base signal that is used in the TDMA signal.
        @param T_slot, slot_occupancy: see package_TDMA_series()
        @param NPSK, NQAM: only one may be > 0, which determines the modulation scheme that's used.
        @return: the signal series.
    """
    def make_signal(t_sample, pattern): # Create the underlying time signal to draw from
        if (NPSK>0):
            samples = PSK_series(t_sample, f_MHz, T_symbol, NPSK)
        else:
            samples = QAM_series(t_sample, f_MHz, T_symbol, NQAM)
        BW = 2/T_symbol # Offset to first null
        samples = scale_signal(samples, 10**(dBm/10.)*1e-3)
        # Assemble into a TDMA frame & suppress out-of-band: typically first sidelobe < -10dB, after suppression < -30dB.
        V_frame = package_TDMA_series(t_sample, [samples]*len(pattern), T_slot, pattern, f_c=f_MHz*MHz,BW=BW,df_stop=BW,stop_dB=20)
        return V_frame
    return pat_generator(make_signal, slot_occupancy, 0, T_slot*len(slot_occupancy))


def generate_IMT(dBm, f_lte_MHz=None, f_ch_MHz=None, ch_slots=["10100001","11111100","01010101"]):
    """
        Constructs International Mobile Telecommunications (i.e. GSM, UMTS) signals.
        These consist of a wide band IMT data channel, plus 200kHz voice channels at each specified channel frequency.
        The time series for each voice channel is composed of frames which consists of 8 timeslots @ 576.92μs per timeslot.
        The 8 timeslots of each frame is ON or OFF as per ch_slots. The frame sequences are repeated to fill time. 
        Voice channels are PSK (symbol duration 3.69231μs), IMT data channels are effectively 256-QAM.
        [ETSI EN 300 910 V8.5.1 (2000-11)]
        @param f_lte_MHz: centre frequency of the IMT data channel [MHz] (default None).
        @param f_ch_MHz: list of centre frequencies of voice channels, matching ch_slots [MHz] (default None).
        @param ch_slots: [slot_occupancy] with slot_occupancy a string consisting of "0" and "1"s, one for each slot.
        @return: the signal series.
    """
    def make_signal(t_sample): # Create the underlying time signal to draw from
        samplerate = 1/np.diff(t_sample).mean()
        V = 0
        if f_lte_MHz: # Wide band LTE data channel.
            # 1.4 - 10 MHz channel widths [https://en.wikipedia.org/wiki/LTE_frequency_bands]
            V = QAM_series(t_sample, f_lte_MHz, 2/3.84e6, 256) # 16 symbols per slot x 16 slots = 256 symbols
            V = scale_signal(V, 10**(dBm/10.)*1e-3)
            # Suppress leakage outside of standard 5 MHz channel width
            V = band_limit(V, samplerate/MHz, [f_lte_MHz-2.5,f_lte_MHz+2.5], [f_lte_MHz-5,f_lte_MHz+5], pass_dB=0.1,stop_dB=20, ftype='cheby1')
        
        if f_ch_MHz: # 200 kHz bandwidth analogue channels
            T_symbol = 3.69231e-6*3 # Time per bit * bits per symbol (NPSK=8)
            T_slot = 156.25*T_symbol
            BW = 2/T_symbol # Offset to first null
            for f_MHz,slots in zip(f_ch_MHz,ch_slots): # Each voice channel and its timeslots arranged into frames
                samples = PSK_series(t_sample, f_MHz, T_symbol, N_symbols=8)
                samples = scale_signal(samples, 10**(dBm/10.)*1e-3)
                # Suppress out-of-band: typically first sidelobe < -10dB, after suppression < -30dB. Spec: 0.5dB ripple over +/-100kHz offset, -30dB at 200kHz offset
                frame = package_TDMA_series(t_sample, [samples]*len(slots), T_slot, slots, f_c=f_MHz*MHz,BW=BW,df_stop=BW,stop_dB=20)
                V += frame
        
        return V
    return sig_generator(make_signal) # TODO: convert to pat_generator with ch_slots - currently ch_slots get re-started on each new generator cycle



############# Unit Testing ####################################################

def assert_isclose(expect, query, message, atol=1e-26, rtol=1e-9):
    expect, query = np.atleast_1d(expect), np.atleast_1d(query)
    if not np.all(np.isclose(query, expect, rtol, atol, equal_nan=True)):
        print("FAILED: %s.\n\tExpected %s, got %s (%s)"%(message,expect,query,str(query-expect)))


def _psd_(f_sample, signal, Navg=1):
    """ Computes the single-sided power spectral density of the signal.
        @return (frequency vector [Hz], power spectral density vector [W/Hz]) """
    fftlen = int(len(signal)/Navg)
    pxx, freqs = mlab.psd(x=signal, NFFT=fftlen, Fs=f_sample, sides="onesided", scale_by_freq=True)
    return freqs, pxx


def _test_PinkNoiseSignal_():
    f_s = 4e9
    t_s = np.arange(2**18)/f_s
    # Test different total power levels
    pass_MHz = [0.2*f_s/MHz,0.3*f_s/MHz]
    stop_MHz = [0.1*f_s/MHz,0.4*f_s/MHz]
    for pdBm in [-100,-60,-30]:
        s = PinkNoiseSignal(t_s, pdBm, pass_MHz, stop_MHz, stop_dB=20)
        psd_f,psd_WHz = _psd_(f_s, s, Navg=256)
        assert_isclose(pdBm, 10*np.log10(np.sum(psd_WHz)*np.diff(psd_f).mean())+30, "PinkNoise sum(PSD) differs from expected power", rtol=0.01)
    
    # Test different pass bands
    pyplot.figure()
    for bw in [0.1,0.2,0.25]:
        pass_MHz = [0.1*f_s/MHz,(0.1+bw)*f_s/MHz]
        stop_MHz = [0.01*f_s/MHz,(0.2+bw)*f_s/MHz]
        s = PinkNoiseSignal(t_s, pdBm, pass_MHz, stop_MHz, stop_dB=20)
        psd_f,psd_WHz = _psd_(f_s, s, Navg=256)
        assert_isclose(pdBm, 10*np.log10(np.sum(psd_WHz)*np.diff(psd_f).mean())+30, "PinkNoise @BW=%g x f_s sum(PSD) differs from expected power"%bw, rtol=0.01)
        pyplot.psd(s, Fs=f_s, NFFT=len(s)//16, scale_by_freq=True, label="BW=%g x f_s"%bw, hold=True)
    pyplot.legend(); pyplot.title("_test_PinkNoiseSignal_")


def _test_GaussianPulsedSignal_():
    f_s = 4e9 # >> 1115MHz, which is used below
    t_s = np.arange(2**19)/f_s
    # Test different total power levels
    for pdBm in [-100,-60,-30]:
        T_width=1e-6; T_separate=10e-6
        s = GaussianPulsedSignal(t_s, pdBm, 1115, T_width, T_separate, pattern="1")
        pulse2carrier_ratio = (0.8*T_width)/T_separate # Fraction of CW power that is contained in a Gaussian pulse that repeats
        psd_f,psd_WHz = _psd_(f_s, s)
        assert_isclose(pdBm, 10*np.log10(1/pulse2carrier_ratio*np.sum(psd_WHz)*np.diff(psd_f).mean())+30,
                       "GaussianPulsedSignal sum(PSD) differs from expected power", rtol=0.01)
    # Test different pulse widths
    pyplot.figure()
    T_separate=10e-6
    for T_width in [4e-6,2e-6,1e-6]:
        s = GaussianPulsedSignal(t_s, pdBm, 1115, T_width, T_separate, pattern="1")
        pulse2carrier_ratio = (0.8*T_width)/T_separate # Fraction of CW power that is contained in a Gaussian pulse that repeats
        psd_f,psd_WHz = _psd_(f_s, s)
        assert_isclose(pdBm, 10*np.log10(1/pulse2carrier_ratio*np.sum(psd_WHz)*np.diff(psd_f).mean())+30,
                       "GaussianPulsedSignal @T_width=%g sum(PSD) differs from expected power"%T_width, rtol=0.01)
        pyplot.plot(t_s/1e-6, s, label="@T_width=%g microsec"%(T_width/1e-6))
    pyplot.ylabel("amplitude"); pyplot.xlabel("time [microsec]")
    pyplot.legend(); pyplot.title("_test_GaussianPulsedSignal_")
    
    # Test different repeat patterns
    pyplot.figure()
    T_width = 1e-6 # Must be short - see CAUTION below
    T_separate = 3*T_width # Repeat pulses truncated at -/+3sigma levels
    for N_off in [0,1,2,3,4,8,10]:
        if (1+N_off > max(t_s)/T_separate /2.): # CAUTION: pulse2carrrier_ratio requires (1+N_off) << max(t_s)/T_separate
            print("INFO: Too short series, skipping test with N_off=%d"%N_off)
            continue
        s = GaussianPulsedSignal(t_s, pdBm, 1115, T_width, T_separate, pattern="1"+("0"*N_off))
        pulse2carrier_ratio = (0.8*T_width)/T_separate * 1/(1.+N_off) # Fraction of CW power that is contained in a Gaussian pulse that repeats
        psd_f,psd_WHz = _psd_(f_s, s)
        assert_isclose(pdBm+10*np.log10(pulse2carrier_ratio), 10*np.log10(np.sum(psd_WHz)*np.diff(psd_f).mean())+30,
                       "GaussianPulsedSignal @1 ON, %d OFF sum(PSD) differs >1 dB from expected power"%N_off, atol=1) # A bit wide, but fails otherwise...
        pyplot.plot(t_s/1e-6, s-N_off*np.max(s), label="@1 ON, %d OFF"%N_off)
    pyplot.ylabel("amplitude"); pyplot.xlabel("time [microsec]")
    pyplot.legend(); pyplot.title("_test_GaussianPulsedSignal_")


import scipy.interpolate
def mcd(t_sample, signal, f_c, BW, downsample=True):
    """
        Determines the in-phase and quadrature phase components of a modulated signal.
        @param t_sample: time series for which to generate the signal [sec].
        @param signal: the series within which the modulated signal is embedded.
        @param f_c: the centre frequency of the modulated signal of interest [Hz].
        @param BW: the bandwidth of the modulated signal of interest [Hz].
        @param downsample: False to leave the I,Q components at the original sampling rate, or True to downsample to 10*BW/2 (default True).
        @return: (in-phase, quadrature phase) signal components
    """
    samplerate = 1/np.diff(t_sample).mean()
    # Band limit
    signal = band_limit(signal, samplerate, [f_c-BW/2.,f_c+BW/2.], [f_c-BW,f_c+BW], pass_dB=0.01, stop_dB=20)
    # Project I & Q components
    s_sin = signal * np.sin(2*np.pi*f_c*t_sample)
    s_cos = signal * np.cos(2*np.pi*f_c*t_sample)
    s_sin = band_limit(s_sin, samplerate, [0,BW/2.], [0,BW], pass_dB=0.01, stop_dB=20)
    s_cos = band_limit(s_cos, samplerate, [0,BW/2.], [0,BW], pass_dB=0.01, stop_dB=20)
    
    if downsample: # Down-sample since the signal is now band-limited
        x = np.arange(0,len(signal), 1)
        x_new = np.arange(0,len(signal), samplerate/(8 * BW/2.)) # Over-sample by factor 8 (>4) to represent trajectory between constellation points
        s_sin = scipy.interpolate.interp1d(x, s_sin, 'cubic')(x_new)
        s_cos = scipy.interpolate.interp1d(x, s_cos, 'cubic')(x_new)
    
    return s_sin, s_cos

def _demo_MOD_series_():
    # Visually confirm bandwidth for constant symbol time
    f_s = 2e9
    t_s = np.arange(2**20)/f_s
    
    f_c = f_s/11.
    T_symbol = (t_s.max()-t_s.min())/4000. # 4000 symbols in the series
    BW = 2/T_symbol # Width to first null
    
    # PSD's
    pyplot.figure()
     
    samples = ASK_series(t_s, f_c/MHz, T_symbol)
    pyplot.subplot(4,1,1)
    pyplot.psd(samples, Fs=f_s, NFFT=len(samples)//16, scale_by_freq=True, label="ASK: BW=%g MHz?"%(BW/1e6))
    pyplot.legend()
     
    pyplot.subplot(4,1,2)
    for NSYMB in [2,8,16]:
        samples = PSK_series(t_s, f_c/MHz, T_symbol, NSYMB)
        pyplot.psd(samples, Fs=f_s, NFFT=len(samples)//16, scale_by_freq=True, label="NPSK=%d: BW=%g MHz?"%(NSYMB,BW/1e6))
    pyplot.legend()
     
    pyplot.subplot(4,1,3)
    for NSYMB in [2,8,16]:
        samples = FSK_series(t_s, f_c/MHz, T_symbol, NSYMB, ff=0.1)
        pyplot.psd(samples, Fs=f_s, NFFT=len(samples)//16, scale_by_freq=True, label="NFSK=%d"%(NSYMB))
    pyplot.legend()
     
    pyplot.subplot(4,1,4)
    for NSYMB in [4,8,16,64]:
        samples = QAM_series(t_s, f_c/MHz, T_symbol, NSYMB)
        pyplot.psd(samples, Fs=f_s, NFFT=len(samples)//16, scale_by_freq=True, label="NQAM=%d: BW=%g MHz?"%(NSYMB,BW/1e6))
    pyplot.legend()
 
    pyplot.title("_demo_MOD_series_")

    # Constellation diagrams
    pyplot.figure()
    
    pyplot.subplot(2,1,1)
    samples = PSK_series(t_s, f_c/MHz, T_symbol, NSYMB)
    pyplot.plot(*mcd(t_s, samples, f_c, 2/T_symbol, downsample=True), '.', label="%d-PSK"%NSYMB)
    pyplot.legend()
     
    pyplot.subplot(2,1,2)
    samples = QAM_series(t_s, f_c/MHz, T_symbol, NSYMB)
    pyplot.plot(*mcd(t_s, samples, f_c, 2/T_symbol, downsample=True), '.', label="%d-QAM"%NSYMB)
    pyplot.legend()
    
    pyplot.title("_demo_MOD_series_")
    

def _test_make_TDMA_slot_():
    f_s = 2e9
    t_s = np.arange(2**16)/f_s
    # Test that power scales as expected
    T_symbol = (t_s.max()-t_s.min())/2000. # Toal number of symbols, need >1000 for 1/N (>=200) to converge to < 10%
    s = ASK_series(t_s, f_s/MHz/4., T_symbol)
    p = signal_power(s)
    for N in [3,4,5]:
        for S in range(N):
            s_slot = _make_TDMA_slot_(t_s, s, T_slot=T_symbol*10, slot_no=S, N_slots=N) # 10 symbols per slot, slots repeat every N*10 symbols
            assert_isclose(p/float(N), signal_power(s_slot), "Slot %d/%d power not scaled as 1/Nslots"%(S,N), rtol=0.1)
    
    # A visual inspection
    pyplot.figure()
    pyplot.psd(s, Fs=f_s, NFFT=len(s)//16, scale_by_freq=True, label="Continuous ASK signal")
    pyplot.psd(s_slot, Fs=f_s, NFFT=len(s)//16, scale_by_freq=True, label="Slot=%d/%d ASK signal"%(S,N), hold=True)
    pyplot.legend(); pyplot.title("_test_make_TDMA_slot_")
    
    pyplot.figure()
    pyplot.plot(t_s/1e-6, s, label="Continuous ASK signal")
    pyplot.plot(t_s/1e-6, s_slot, label="Slot=%d/%d ASK signal"%(S,N))
    pyplot.ylabel("amplitude"); pyplot.xlabel("time [microsec]")
    pyplot.legend(); pyplot.title("_test_make_TDMA_slot_")
    

def show(t_sample, S, title, f_resolution=1e6, nyqzone=1, minmax=True, overplot_fig=None):
    """
        Generates three figures to display a waterfall map (time vs frequency), time series amplitude and also a Power Spectral Density plot.
        The Blackman-Harris windowing function is applied in the FFT.
        @param t_sample: series of sampling instants [s].
        @param S: the (amplitude) signal series sampled at time instants t_sample.
        @param nyqzone: Nyquist zone in which the signal has been sampled (default 1).
        @param minmax: True to add min & max hold to PSD in faint grey lines (default True).
        @param overplot_fig: list of two figure numbers to overplot time & PSD onto (default None).
    """
    S = S*(1e3)**.5 # Convert from W to mW scale
    f_s = 1/np.diff(t_sample).mean()
    f_c = (nyqzone-1)*f_s/2
    if (nyqzone %2 == 0): # Flip sideband for display purposes by mixing with f_s/2
        S[::2] *= -1
    
    fftlen = (int(f_s/f_resolution) //2)*2  # Ensure it's at least an even number, or else FFT is extremely slow
    win = scipy.signal.windows.blackmanharris(fftlen, sym=False) # sym=False for a DFT-even window as recommended in [Harris 1976/78]
    win /= np.sum(np.abs(win)**2)**.5 # Normalized for processing gain of window
    
    # A spectrogram
    N_spectra = int(len(S) / fftlen)
    V = np.reshape(S[:N_spectra*fftlen],(-1,fftlen)) # Break into time chunks
    FV = np.fft.fft(V*win, axis=1, n=fftlen) * (1/f_s)**.5 # Normalize numpy FFT to "/sqrt(Hz)"
    FV = 2*np.abs(FV[:,:fftlen//2])**2
    FV[0] /= 2. # Half the DC bin
    FV = 10*np.log10(FV) # To "Power spectral density [dBm/Hz]"
    freq = f_c + np.fft.fftfreq(fftlen, 1.0/f_s)[:fftlen//2]
    time = np.arange(N_spectra)*fftlen/f_s + t_sample[0]
    freq /= 1e6 # MHz
    time /= 1e-6 # microsec
    pyplot.figure()
    pyplot.imshow(FV, aspect='auto', extent=(freq[0],freq[-1], time[-1],time[0]),
                  vmin=np.percentile(FV,66), vmax=np.percentile(FV,99), cmap='hot', interpolation=None)
    pyplot.xlabel('Frequency [MHz]'); pyplot.ylabel('Time [$\mu$s]')
    pyplot.title(title)

    # A time series
    STEP = int(1 + len(S)//1e5) # This limits figures to no more than 1e5 points
    pyplot.figure(overplot_fig[0] if overplot_fig else None)
    pyplot.plot(t_sample[::STEP]/1e-6, S[::STEP], '.-', alpha=0.4)
    pyplot.xlabel('time [microsec]'); pyplot.ylabel('Amplitude [V/sqrt(Ohm)]'); pyplot.grid(True)
    pyplot.title(title)

    # A PSD
    pyplot.figure(overplot_fig[1] if overplot_fig else None)
    if minmax:
        pyplot.plot(freq, np.max(FV,axis=0), 'k', freq, np.min(FV,axis=0), 'k', alpha=0.1)
    pyplot.plot(freq, np.mean(FV,axis=0))
    pyplot.xlabel('Frequency [MHz]'); pyplot.ylabel('PSD [dBm/Hz]'); pyplot.grid(True)
    pyplot.title(title)

    
def _demo_RFIUSECASE_(f_s=16e9, n_pts=2**20):
    """
        Demo of some RFI Use Case, covering 100 - 4000 MHz.
        Generates a PSD and a spectrogram
    """
    t_sample = np.arange(n_pts)/f_s
    
    # Some signals
    S = generate_DME(-55, 1050, pattern="010100")+generate_DME(-55, 1115, pattern="100100")
    S = S(t_sample) + WhiteNoiseSignal(t_sample, Teq=0.1) # Add some receiver noise
    show(t_sample, S, "DME")
    
    # More signals
    S = generate_TDMA(-100, 235, 2/20e3, 1, "1") # DAB, always on
    S += generate_TDMA(-100, 680, 2/8e6, 1, "1", NPSK=8) # DTTV, always on
    S += generate_IMT(-80, 930, [945,950],["10","11"]) # IMT-900 with two voice channels, may be on simultaneously
    S += generate_SSR(-55, T_repeat=55e-6)
    S += generate_DME(-55, 1115, pattern="100100") +\
         generate_DME(-55, 1050, pattern="100010").delay(6e-6) # For demo peaks DON'T coincide -- otherwise typically leads to clipping, which should just result in time-domain excision.
    S += generate_TDMA(-85, 1250, 2/10e6, 1, "1") # GNSS, always on
    S += np.sum([generate_TDMA(-70, 1545+i/2., 2/10.5e3, 1, "1", NPSK=4) for i in range(10)], axis=0) # INMARSAT Aero Sat-to-Plane
    S += np.sum([generate_TDMA(-105, 1620+i, 2/31.5e3, 128/31.5e3, "1010", NPSK=4) for i in range(10)], axis=0) # IRIDIUM
    S += generate_IMT(-90, 1950) # IMT-1900, always on
    S += generate_TDMA(-60, 3150, 2/500e3, 1, "1100") # RNAV
    S = S(t_sample) + WhiteNoiseSignal(t_sample, Teq=30) # Add some receiver noise
    show(t_sample, S, "_demo_RFIUSECASE_", f_resolution=.25e6) # Resolves INMARSAT & IRIDIUM combs at expense of time resolution
    

    
if __name__ == "__main__":
    _test_PinkNoiseSignal_()
    _test_GaussianPulsedSignal_()
    _demo_MOD_series_()
    _test_make_TDMA_slot_()
    print("ALL TESTS COMPLETED")
    _demo_RFIUSECASE_()
