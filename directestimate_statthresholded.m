%   adapted from the code written for Kramer et al. (2008) 
%   by Tolga Esat Ozkurt (2011)
%   to implement direct PAC estimate (?zkurt and Schnitzler, 2011)
%   and raw moduation index (Canolty et el., 2006)  
%   statistically thresholded as suggested by (?zkurt, 2012) 
%
%   INPUTS:
%    d      = The unfiltered data.
%    srate  = The sampling rate (e.g., 1000 Hz)
%    r_value = 1 - p_value ~ 0.99
%
%   OUTPUS:
%    
%    mod2d  = The two-dim modulation index.
%    flow   = The frequency axis for the lower phase frequencies.
%    fhigh  = The frequency axis for the higher amplitude frequencies.
%
%   For example use the command: [DE MI CC aa bb] = directestimate_statthresholded(data, fs,2:2:60,100:10:400,4,10,r,1);
%   You may take r = 1 - (0.01 / (30*30)); for the example above in order to compensate for multiple comparison (bonferroni correction) or take r = 0.99
%
%   Please cite these two references in case you use this routine.
%
%   T. E. ?zkurt, A. Schnitzler, ?A critical note on the definition of phase-amplitude cross-frequency coupling?, Journal of Neuroscience Methods, vol. 201, no. 2, pp. 438-443, 2011.
%   T. E. ?zkurt, "Statistically reliable and fast direct estimation of phase-amplitude cross-frequency coupling", IEEE Transactions on Biomedical Engineering, vol. 59. no. 7, pp. 1943-1950, 2012.


 
function [mod2dDirectEst_thr, mod2dDirectEst, modConf, flow, fhigh] = directestimate_statthresholded(d, srate,flow1,fhigh1,bwlow,bwhigh,r_value,plotfigure)
 
  %The phase frequencies.
  flow2 = flow1+bwlow;                               
  
  %The amp frequencies.
  fhigh2 = fhigh1+bwhigh;                           
  
  mod2dRawMI = zeros(length(flow1), length(fhigh1));
  mod2dDirectEst = zeros(length(flow1), length(fhigh1));
  modConf = zeros(length(flow1), length(fhigh1));
  
  for i=1:length(flow1)
      theta=eegfilt(d,srate,flow1(i),flow2(i));  %Compute the low freq signal.
      theta=theta(srate:length(d)-srate-1);           %Drop the first and last second.
      phase = angle(hilbert(theta));                  %Compute the low freq phase.
      ['Loops remaining = ', num2str(length(flow1)-i)]
      for j=1:length(fhigh1)
        gamma=eegfilt(d,srate,fhigh1(j),fhigh2(j));%Compute the high freq signal.
        gamma=gamma(srate:length(d)-srate-1);         %Drop the first and last second.
        amp = abs(hilbert(gamma));                 %Compute the high freq amplitude.
        
        amp = (amp - mean(amp))./std(amp); % normalization
       
        
        %Compute the modulation index.
        [m_norm1 m_norm2 m_norm3] = modulation_index_nostats(amp, phase);
        mod2dRawMI(i,j) = m_norm1;
        mod2dDirectEst(i,j) = m_norm2;
        modConf(i,j) = m_norm3;
      end
  end
 
  %% Plot the two-dimensional PAC portraits
  
  flow = (flow1 + flow2) / 2;
  fhigh = (fhigh1 + fhigh2) / 2;


  low_limit = 2*((erfinv(r_value)*sqrt(length(d)))^2);
  mod2dDirectEst_thr = mod2dDirectEst;
  mod2dDirectEst_thr(modConf < low_limit) = 0;
  
  if plotfigure 
      
      

      figure
      imagesc(flow, fhigh, mod2dDirectEst_thr');  colorbar;
      axis xy
      set(gca, 'FontSize', 12);
      xlabel('Phase Frequency [Hz]');  ylabel('Amplitude Frequency [Hz]');
      
      figure
      imagesc(flow, fhigh, modConf');  colorbar;
      axis xy
      set(gca, 'FontSize', 12);
      xlabel('Phase Frequency [Hz]');  ylabel('Amplitude Frequency [Hz]');
      
  end


  
  %% end of the function.
  
  
 
  function [m_raw1 m_raw2 m_raw3] = modulation_index_nostats(a, p)
 
  N = length(a);
  z = a.*exp(1i*p);
  
  m_raw1 = abs(mean(z));   %Compute the mean length 
  m_raw2 = (1./sqrt(N)) * abs(sum(z)) / sqrt(sum(a.*a)); % compute the direct estimate
  m_raw3 = abs(sum(z)) ^ 2;