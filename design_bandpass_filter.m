% Designs a bandpass filter.m
%
% INPUTS
%
% stopband_freq1     (Hz>0). Frequencies at and below this value are to be attenuated
% passband_freq1     (Hz>0). Lower limit of frequency to pass through filter with little or no attenuation
% passband_freq2     (Hz>0). Upper "                                      
% stopband_freq2     (Hz>0). Frequencies at and above this value are to be attenuated
% stopband_atten1    (dB>0). Attenuation in the stop band at frequencies <=stopband_freq1
% passband_ripple    (dB>0). Maximum ripple in pass band. reasonable value is 1.
% stopband_atten2    (dB>0). Attenuation in stop band at frequencies >=stopband_freq2
% design_method      Ascii string. Could use 'ellip'. 
%                    The following is from matlab documentation for R2016b:
%                    'DesignMethod'  ? Design method'butter' | 'cheby1' | 'cheby2' | 'cls' | 'ellip'  | 'equiripple' | 'freqsamp' | 'kaiserwin' | 'ls'  | 'maxflat' | 'window'
%
%                   Design method, specified as the comma-separated pair consisting of 'DesignMethod' and a string. The choice of design method depends on the set of frequency and magnitude constraints that you specify.
% 
%                     *
% 
%                       'butter' designs a Butterworth IIR filter. Butterworth filters have a smooth monotonic frequency response that is maximally flat in the passband. They sacrifice rolloff steepness for flatness.
%                     *
% 
%                       'cheby1' designs a Chebyshev type I IIR filter. Chebyshev type I filters have a frequency response that is equiripple in the passband and maximally flat in the stopband. Their passband ripple increases with increasing rolloff steepness.
%                     *
% 
%                       'cheby2' designs a Chebyshev type II IIR filter. Chebyshev type II filters have a frequency response that is maximally flat in the passband and equiripple in the stopband.
%                     *
% 
%                       'cls' designs an FIR filter using constrained least squares. The method minimizes the discrepancy between a specified arbitrary piecewise-linear function and the filter's magnitude response. At the same time, it lets you set constraints on the passband ripple and stopband attenuation.
%                     *
% 
%                       'ellip' designs an elliptic IIR filter. Elliptic filters have a frequency response that is equiripple in both passband and stopband.
%                     *
% 
%                       'equiripple' designs an equiripple FIR filter using the Parks-McClellan algorithm. Equiripple filters have a frequency response that minimizes the maximum ripple magnitude over all bands.
%                     *
% 
%                       'freqsamp' designs an FIR filter of arbitrary magnitude response by sampling the frequency response uniformly and taking the inverse Fourier transform.
%                     *
% 
%                       'kaiserwin' designs an FIR filter using the Kaiser window method. The method truncates the impulse response of an ideal filter and uses a Kaiser window to attenuate the resulting truncation oscillations.
%                     *
% 
%                       'ls' designs an FIR filter using least squares. The method minimizes the discrepancy between a specified arbitrary piecewise-linear function and the filter's magnitude response.
%                     *
% 
%                       'maxflat' designs a maximally flat FIR filter. These filters have a smooth monotonic frequency response that is maximally flat in the passband.
%                     *
% 
%                       'window' designs an FIR filter using the windowing method. Window
% fs                 Sample frequency (Hz)
% fign               Matlab figure window to plot filter design. If empty, no plot made
%
% OUTPUTS
%

function [d]=design_bandpass_filter(stopband_freq1,passband_freq1,passband_freq2,stopband_freq2,...
    stopband_atten1,passband_ripple,stopband_atten2,design_method,fs,fign)

% save jlsee

d = designfilt('bandpassiir', ...       % Response type
    'StopbandFrequency1',stopband_freq1, ...    % Frequency constraints
    'PassbandFrequency1',passband_freq1, ...
    'PassbandFrequency2',passband_freq2, ...
    'StopbandFrequency2',stopband_freq2, ...
    'StopbandAttenuation1',stopband_atten1, ...   % Magnitude constraints
    'PassbandRipple',passband_ripple, ...
    'StopbandAttenuation2',stopband_atten2, ...
    'DesignMethod',design_method, ...      % Design method
    'MatchExactly','passband', ...   % Design method options
    'SampleRate',fs);               % Sample rate



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the filter design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(fign))
    % set # frequencies to plot
    nfreq=1000;
    % Get highest plotted frequency (Hz)
    fmax=fs/2;
    % Get frequencies at which filter design is plotted (Hz).
    f=linspace(0,fmax,nfreq);
    
    h=freqz(d,f,fs);
    
    % Get magnitude and phase
    a=abs(h);
    p=atan2(imag(h),real(h));
    % Get attenuation in dB
    db=20*log10(a);
    
    % Set min/max values of attenuation (dB) for plot
    max_db=2;
    min_db=-70;
 
    figure(fign);
    subplot(2,1,1);
    plot(f,db,'b-');
    V=axis;
    V([3 4])=[min_db max_db];
    axis(V);
    grid on
    xlabel('Hz');
    ylabel('dB');
    title('Filter Design');
    subplot(2,1,2);
    plot(f,p,'r-');
    grid on
    xlabel('Hz');
    ylabel('Radians');
end % if(~isempty(fign))