close all; clear;

tau = 2.^(4:7);
dV = 10.^(-3.5:0.5:-2);

K = [18246753246.753246, 15056818181.818182, 5754870129.870129, 17467532467.532467; ...
     7288961038.961037, 5438311688.311686, 8262987012.987013, 6095779220.779219; ...
     5608766233.766233, 5024350649.350649, 5024350649.350649, 3928571428.5714264; ...
     4172077922.077923, 4488636363.636364, 3831168831.168831, 4318181818.181818];  % Nv x Nt
 
d_K = [[1, 1, 1, 1] * 30; ...
       [1, 1, 1, 1] * 10; ...
       [1, 1, 1, 1] * 3; ...
       [1, 1, 1, 1] * 1] * 1.25e10 / 92;
   
N_dV = length(dV);
N_tau = length(tau);
 
getFig('$\tau$ (ns)', '$dV/V$', '$K$ (Pa)', 'log', 'log', 'log');
surf(tau, dV, K, 'EdgeColor', 'interp', 'FaceColor', 'interp');

getFig('$\tau_T$ (ns)', '$K$ (GPa)', '$K(\tau_T, dV/V)$', 'log', 'log');
for i = 1:N_dV
    errorbar(tau, K(i, :) / 1e9, d_K(i, :) / 1e9, 'o', ...
        'DisplayName', ['dV/V = 1e' num2str(log10(dV(i)))], 'LineWidth', 1);
end
xlim([9.9, max(tau, [], 'all') * 1.1]);
ylim([min(K, [], 'all') * 0.9, max(K, [], 'all') * 1.1] / 1e9);

getFig('$dV/V$', '$K$ (GPa)', '$K(\tau_T, dV/V)$', 'log', 'log');
for i = 1:N_tau
    errorbar(dV, K(:, i) / 1e9, d_K(:, i) / 1e9, 'o', ...
        'DisplayName', ['tau = ' num2str(tau(i))], 'LineWidth', 1);
end
xlim([min(dV, [], 'all') * 0.9, max(dV, [], 'all') * 1.1]);
ylim([min(K, [], 'all') * 0.9, max(K, [], 'all') * 1.1] / 1e9);
 