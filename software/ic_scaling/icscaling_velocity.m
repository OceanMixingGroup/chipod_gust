function [ eps, var_eps, eps_f] = icscaling_velocity(spec_time, spec_f, spec, vel_time, vel, f_range)
%% [ eps, var_eps, eps_f] = icscaling_velocity(spec_time, spec_f, spec, vel_time, vel, f_range)
%  
%     This function calculates epsilon based on a inertail convective scaling (ic) 
%     of a velocity spectrum in a given freuquency range
%
%     INPUT
%        spec_time   :  time vector of spectrogram
%        spec_f      :  frrequency vector of spectrogram
%        spec        :  spectrogram( f, spec_time )
%        vel_time    :  time vector of velocity
%        vel         :  velocity( vel_time )
%        f_range     :  frequency range to perform the ic-scaling
%
%     OUTPUT
%        eps         :  epsilon( spec_time ) 
%        var_eps     :  var(log10( eps_f( f_range) ))
%        eps_f       :  normalized spectrogram unit epsilon eps_f( f, spec_time )
%
%
%   created by: 
%        Johannes Becherer
%        Tue Oct 10 13:21:06 PDT 2017


%_____________________check dimensions______________________
if ~(size(spec,1) == length(spec_f))
   if ( size(spec,2) == length(spec_f) ) & ( size(spec,1) == length(spec_time) )
      spec = spec';
   else
      error('The size of spec does not matach the size of spec_f and spec_time')
   end
end



%_____________________get vel( spec_time )______________________

 [vel_spec]   = clever_interp( vel_time, vel, spec_time );


%_____________________normalize spectrum______________________
  % C = 11/15;
  C = .6;
  
  eps_f =  (2*pi/C^1.5.*spec.^1.5.*spec_f.^2.5)./abs(vel_spec);

%_____________________average eps time series in f_range______________________
		ii_fband = find( spec_f>=f_range(1) &  spec_f<=f_range(2) );
		eps      = mean( eps_f( ii_fband,: ));
      var_eps  = var( log10(eps_f( ii_fband,: )));


