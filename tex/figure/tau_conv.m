close all; clear;

tau_20 = 2.^(2:9);
tau_40 = 2.^(3:2:9);
ids_40 = 2:4;
ids_20 = 4:8;

K_my_20 = [1428227011.3276384, 1999297828.012307, 2835870424.072775, 3741000507.362456, 3691978215.1999416, 4650613874.690208, 4382595767.6775, 5668130528.288169];
K_gmx_20 = [1047548803.8262933; 1203164355.4879146; 1118967446.186921; 1027025695.0268532; 686843623.7888445; 503773007.2384608; 260497180.7873671; 92492176.13782674];
K_my_40 = [2003652026.6738977, 3587152248.178136, 3884408276.698165, 3819830967.5079246];
K_gmx_40 = [1513239263.0459206, 1564820691.4133577, 934654343.1573989, 197501832.05088243];

d_K_my_20 = K_my_20 / sqrt(20 / 0.005);
d_K_gmx_20 = K_gmx_20 / sqrt(20 / 0.005);
d_K_my_40 = K_my_40 / sqrt(40 / 0.005);
d_K_gmx_40 = K_gmx_40 / sqrt(40 / 0.005);

K0_40 = mean(K_my_40(ids_40));
Kfit_40 = polyfit(log(tau_40(ids_40)), log(K_my_40(ids_40)), 1);
Kfit_20 = polyfit(log(tau_20(ids_20)), log(K_my_20(ids_20)), 1);

fig_40 = getFig('$\tau_P$ (ns)', 'K (GPa)', '$K(\tau_P, \tau_T = \tau_P / 4)$, $t_{avg} = 40$ ns', 'log', 'log');
plot(fig_40.ax, tau_40, K_my_40, 'o', 'DisplayName', 'my analysis', 'Color', getMyColor(1));
plot(fig_40.ax, tau_40, K_gmx_40, 'o', 'DisplayName', 'gromacs');
plot(fig_40.ax, [min(tau_40), max(tau_40)], [1, 1] * K0_40, '--', 'DisplayName', ['$K_{stab} = ' num2str(round(K0_40 / 1e9, 2)) '$ GPa'], 'Color', getMyColor(1));
y_lims_40 = get(fig_40.ax, 'YLim');
ylim(fig_40.ax, [y_lims_40(1), 5e9]);

fig_20 = getFig('$\tau_P$ (ns)', 'K (GPa)', '$K(\tau_P, \tau_T = \tau_P / 4)$, $t_{avg} = 20$ ns', 'log', 'log');
plot(fig_20.ax, tau_20, K_my_20, 'o', 'DisplayName', 'my analysis', 'Color', getMyColor(1));
plot(fig_20.ax, tau_20, K_gmx_20, 'o', 'DisplayName', 'gromacs');
plot(fig_20.ax, tau_20, exp(polyval(Kfit_20, log(tau_20))), '--', 'DisplayName', '$\neq$ const', 'Color', getMyColor(1));
