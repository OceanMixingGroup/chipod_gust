function [W] = pitot_choose_W(rdat)
% W = pitot_choose_W(rdat)
%
% This code carefully chooses whether to use W or WP as the input for the
% pitot. In newer chipods and gusTs, the correct data is almost always
% saved in W, however, there are some older chipods (2015ish) where W and WP
% were switched and the correct data was saved in WP. (In VERY old chipods from
% the Dynamo era, W and WP were saved at W2 and W3.)
%
% The goal of this code is to remain backwards compatable with old chipods
% while reducing the instances where WP is incorrectly chosen as the pitot
% voltage instead of W.
%
% This code fixes three issues:
% 1. In new pitot circuit boards (created by Pavan in 2018), WP is saved as
% near-zero volts, so the old method of choosing WP vs W based on which had
% a mean that was further from 2.02 doesn't work with these newer boards.
%
% 2. In cases where the good pitot voltage is close to 2 volts (usually 
% occurs when velocities are > 1 m/s), the code may incorrectly choose WP 
% rather than W. I've added an error message, so that the user can see when 
% this is happening and manually hard-wire the code to choose W instead of WP.
%
% 3. I have put the algorithm that chooses WP vs W into this function 
% rather than having it hard-line coded into the following functions: 
%           chipod_gust/software/pitot/pitot_avg_raw_data.m
%           chipod_gust/software/pitot/chi_calibrate_chipod_pitot.m
%           chipod_gust/software/main_proc/chi_calibrate_chipod.m
%           chipod_gust/software/pitot/chi_calibrate_gust_pitot.m
%           chipod_gust/software/main_proc/chi_calibrate_gust.m
%
% written by Sally Warner, June 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you are processing a chipod that was deployed since 2015 and you
% received the warning message about the code choosing WP instead of W,
% please uncomment line 40 (W  = rdat.W;) and comment out everything below.

% W  = rdat.W;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(rdat, 'W')
    
    % Case using pitot circuit boards newer than 2018 where WP is near
    % zero. Always want to use W, not WP, as the pitot voltage.
    if nanmean(rdat.WP) < 10^-3
        W  = rdat.W;
        
    % case using older pitot circuit boards where the average of WP is
    % close to 2.02. Code leaves room to choose WP in the SUPER
    % RARE cases where the wires to WP and W were switched (circa 2015).
    else
    
        dV1 = abs(nanmean(rdat.W)-2.02);
        dV2 = abs(nanmean(rdat.WP)-2.02);
        if dV1>dV2
            W  = rdat.W;
        else
            W  = rdat.WP;
            disp({'*** WARNING *** pitot_choose_W.m is choosing WP instead';...
                'of W as the pitot voltage. This is only the correct choice in ';...
                'VERY RARE INSANCES (mostly around 2015). Likely, you should';...
                'hard-wire the code to always choose W as the pitot voltage in:';...
                'chipod_gust/software/pitot/pitot_choose_W.m'})
        end
    end
    
% Very very old case (circa Dynamo) where W and WP were saved as W2 and W3.  
else
	dV1 = abs(nanmean(rdat.W2)-2.02);
	dV2 = abs(nanmean(rdat.W3)-2.02);
	if dV1>dV2
        W  = rdat.W2;
	else
        W  = rdat.W3;
	end
end


