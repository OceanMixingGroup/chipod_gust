function [do] = ChooseEstimates(do, commands)
% Choose estimates to process
%       [do] = ChooseEstimates(do, commands)
% 1. With no arguments, returns defaults
% 2. Available commands:
%         'no_vc': no viscous-convective estimates
%         'no_ic': no inertial-convective estimates
%         'no_pitot': no pitot estimates
%         'none': reset all to 0 [EXAMPLE: do = ChooseEstimates([], 'none')]
% 3. commands can be chained: 'no_ic; no_pitot'

    if ~exist('commands', 'var'), commands = ''; end
    if ~exist('do', 'var') | isempty(do)
        do = [];
        % turn off vc estimates
        do.chi_mi11     = 1;
        do.chi_mi22     = 1;
        do.chi_mm1      = 1;
        do.chi_mm2      = 1;
        do.chi_pi11     = 1;
        do.chi_pi22     = 1;
        do.chi_pm1      = 1;
        do.chi_pm2      = 1;

        % turn off ic estimates
        do.chi_mi11_ic  = 1;
        do.chi_mi22_ic  = 1;
        do.chi_mm1_ic   = 1;
        do.chi_mm2_ic   = 1;
        do.chi_pi11_ic  = 1;
        do.chi_pi22_ic  = 1;
        do.chi_pm1_ic   = 1;
        do.chi_pm2_ic   = 1;

        % turn off cases where 1&2 was used for dT/dz
        do.chi_mi112	= 0;
        do.chi_mi212	= 0;
        do.chi_pi112	= 0;
        do.chi_pi212	= 0;
        do.chi_mi112_ic	= 0;
        do.chi_mi212_ic     = 0;
        do.chi_pi112_ic	= 0;
        do.chi_pi212_ic     = 0;
    end

    if ~isempty(strfind(commands, 'no_vc'))
        % turn off vc estimates
        do.chi_mi11     = 0;
        do.chi_mi22     = 0;
        do.chi_mm1      = 0;
        do.chi_mm2      = 0;
        do.chi_pi11     = 0;
        do.chi_pi22     = 0;
        do.chi_pm1      = 0;
        do.chi_pm2      = 0;
    end

    if ~isempty(strfind(commands, 'no_ic'))
        % turn off ic estimates
        do.chi_mi11_ic  = 0;
        do.chi_mi22_ic  = 0;
        do.chi_mm1_ic   = 0;
        do.chi_mm2_ic   = 0;
        do.chi_pi11_ic  = 0;
        do.chi_pi22_ic  = 0;
        do.chi_pm1_ic   = 0;
        do.chi_pm2_ic   = 0;
    end

    if ~isempty(strfind(commands, 'no_pitot'))
        do.chi_pm1 = 0;
        do.chi_pm2 = 0;
        do.chi_pm1_ic = 0;
        do.chi_pm2_ic = 0;

        do.chi_pi11 = 0;
        do.chi_pi22 = 0;
        do.chi_pi11_ic = 0;
        do.chi_pi22_ic = 0;
    end

    if ~isempty(strfind(commands, 'no_pm'))
        do.chi_pm1 = 0;
        do.chi_pm2 = 0;
        do.chi_pm1_ic = 0;
        do.chi_pm2_ic = 0;
    end


    if ~isempty(strfind(commands, 'no_pi'))
        do.chi_pi11 = 0;
        do.chi_pi22 = 0;
        do.chi_pi11_ic = 0;
        do.chi_pi22_ic = 0;
    end


    if ~isempty(strfind(commands, 'no_mm'))
        do.chi_mm1 = 0;
        do.chi_mm2 = 0;
        do.chi_mm1_ic = 0;
        do.chi_mm2_ic = 0;
    end


    if ~isempty(strfind(commands, 'no_mi'))
        do.chi_mi11 = 0;
        do.chi_mi22 = 0;
        do.chi_mi11_ic = 0;
        do.chi_mi22_ic = 0;
    end

    if ~isempty(strfind(commands, 'none'))
        % turn off vc estimates
        do.chi_mi11     = 0;
        do.chi_mi22     = 0;
        do.chi_mm1      = 0;
        do.chi_mm2      = 0;
        do.chi_pi11     = 0;
        do.chi_pi22     = 0;
        do.chi_pm1      = 0;
        do.chi_pm2      = 0;

        % turn off ic estimates
        do.chi_mi11_ic  = 0;
        do.chi_mi22_ic  = 0;
        do.chi_mm1_ic   = 0;
        do.chi_mm2_ic   = 0;
        do.chi_pi11_ic  = 0;
        do.chi_pi22_ic  = 0;
        do.chi_pm1_ic   = 0;
        do.chi_pm2_ic   = 0;

        % turn off cases where 1&2 was used for dT/dz
        do.chi_mi112	= 0;
        do.chi_mi212	= 0;
        do.chi_pi112	= 0;
        do.chi_pi212	= 0;
        do.chi_mi112_ic	= 0;
        do.chi_mi212_ic     = 0;
        do.chi_pi112_ic	= 0;
        do.chi_pi212_ic     = 0;
    end
