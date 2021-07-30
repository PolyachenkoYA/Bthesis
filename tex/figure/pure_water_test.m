clear; close all;

t = 2.^(-2:1);
units = 1e3;

K_dV = [2.0303030303030307, 1.987878787878788, 1.9737373737373738, 1.9414141414141417] * units;
K_f = [2.0828282828282827, 2.4101010101010103, 1.9232323232323234, 2.0060606060606063] * units;
d_K_dV = ([10, 4.5, 5, 5] / 88.5 / 2) * units;
d_K_f = ([24, 49, 29.5, 59] / 88.5 / 2) * units;
K0 = 2.098989898989899 * units;

getFig('$t_{stab}$ (ns)', '$K$ (MPa)', 'pure tip4p/2005 water; $T = 30$ $C^{\circ}$; $K(t_{stab})$');
errorbar(t, K_f, d_K_f, 'o', ...
         'DisplayName', '$\langle \Delta V^2 \rangle$ method', 'LineWidth', 1.5);
errorbar(t, K_dV, d_K_dV, 'o', ...
        'DisplayName', '$dP/dV$ method', 'LineWidth', 1.5);
plot([min(t), max(t)], [1, 1] * K0, '--', ...
        'DisplayName', 'experiment value', 'LineWidth', 1.5, 'Color', 'black');
xlim([0, 2.1]);
