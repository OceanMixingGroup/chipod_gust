% compare winters & d'asaro estimate with others
load ../proc/Turb.mat
t0 = find_approx(Turb.mm1.time, chi.wda.time(1));
t1 = find_approx(Turb.mm1.time, chi.wda.time(end));

sgn = sign(interp1(chi.time, chi.dTdz, chi.wda.time));
sgn(sgn == 0) = 1;

load ../input/dTdz_i.mat
tz0 = find_approx(Tz_i.time, chi.wda.time(1));
tz1 = find_approx(Tz_i.time, chi.wda.time(end));

figure('Position', [1 5 1280 1340]);
ax(1) = subplot(311); hold on;
plot(chi.wda.time, chi.wda.Jqi, 'color', [0.8, 0.8, 0.8]);
hwd = plot(chi.wda.time, chi.wda.Jq .* sgn);
plot(Turb.mm1.time, Turb.mm1.Jq, 'k', 'linewidth', 2)
ylabel('Jq')
legend('internal', 'W&DA', 'mooring')
xlim(chi.wda.time([1 1440]))
ylim([prctile(chi.wda.Jq .* sgn, 1), prctile(chi.wda.Jq .* sgn, 99.5)]);
datetick('keeplimits')

ax(2) = subplot(312);
semilogy(chi.wda.time, chi.wda.Kt, 'color', hwd.Color); hold on;
semilogy(Turb.mm1.time, Turb.mm1.Kt, 'k', 'linewidth', 2)
ylabel('Kt')
ylim([1e-9, 1])
xlim(chi.wda.time([1 1440]))
datetick('keeplimits')

ax(3) = subplot(313);
plot(Tz_i.time, Tz_i.Tz1, 'color', [0.8, 0.8, 0.8]); hold on;
plot(chi.wda.time, chi.wda.dTdz .* sgn);
plot(Turb.mm1.time, Turb.mm1.dTdz, 'k', 'linewidth', 2)
ylabel('dTdz')
linkaxes(ax, 'x')
datetick('keeplimits')
xlim(chi.wda.time([1 1440]))
plot(xlim, [0, 0], 'color', [1 1 1]*0.75);
