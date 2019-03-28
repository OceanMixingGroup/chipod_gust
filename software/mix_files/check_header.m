% will create header.mat file if necessary
% if header.mat exists, it will read it
head = chi_get_calibration_coefs(basedir);

% will create header.mat file if necessary
% if header_p.mat exists, it will read it
W  = chi_get_calibration_coefs_pitot(basedir);

disp(' ')
disp('----- temp T1; T2')
[head.coef.T1'; head.coef.T2']

disp(' ')
disp('----- differentiator T1P; T2P')
[head.coef.T1P'; head.coef.T2P']

disp(' ')
disp('----- Accelerometer: AX; AY: AZ')
[head.coef.AX'; head.coef.AY'; head.coef.AZ']

disp(' ')
disp('----- Compass')
disp('First element might have been corrected for declination.')
head.coef.CMP'

disp(' ')
disp('----- Pitot')
W
