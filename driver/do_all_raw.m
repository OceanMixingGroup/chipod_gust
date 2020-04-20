% do_all_raw.m
%
% I want to be able to run:
%   - do_temp_proc
%   - do_raw_pitot
%   - do_dTdz_i_proc
% in sequence on the matlab servers without having to start each one when
% the previous code finishes.
%
% I couldn't figure out a way to do them in sequence from linux. I found
% a help page that says to use the "&" command to run multiple commands one
% at a time, but that doesn't seem to work. My guess is that there's
% something about running matlab in the background that makes it think a
% process is done when it isn't.
%
% The way to get around this is to just run all three from this simple 
% matlab script rather than from the command line in linux.
%
% April 2020 SJW

clear

try
    do_raw_pitot
catch
    disp('do_raw_pitot crashed')
end

try
    do_temp_proc
catch
    disp('do_temp_proc crashed')
end

try
    do_dTdz_i_proc
catch
    disp('do_dTdz_i_proc crashed')
end