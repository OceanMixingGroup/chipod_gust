% Easy way to open compass file

[filename, pathname] = uigetfile('*.*')

[data, head] = raw_load_chipod([pathname '/' filename])

figure;
plot(data.CMP/10)
ylim([-20, 380])
axis square