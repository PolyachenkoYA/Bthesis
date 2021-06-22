close all; clear;

TC2K = 273.15;

K_Gorelov = [275.15, 1331.9148936170213; ...
277.15, 1031.9148936170213; ...
281.15, 306.3829787234042; ...
291.15, 146.80851063829778; ...

294.7004608294931, 174.71698113207548; ...
298.0184331797235, 159.99999999999999; ...
299.95391705069125, 151.69811320754718; ...
302.7188940092166, 139.62264150943395; ...
306.1290322580645, 151.32075471698114; ...
309.7235023041475, 107.9245283018868; ...
319.5852534562212, 95.47169811320755; ...
325.852534562212, 35.47169811320754];


getFig('$T$ $C^{\circ}$', '$K$ MPa', 'K(T)', '', 'log');
plot(K_Gorelov(:, 1) - TC2K, K_Gorelov(:, 2), 'o');
