function K1 = get_K_interp(T1, verbose)
    TC2K = 273.15;

    Gorelov_data = [275.15, 1331.9148936170213; ...
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

    T_b = 20;
    T = Gorelov_data(:, 1) - TC2K;
    K = Gorelov_data(:, 2);

    if(verbose)
        fig = getFig('$T$ ($C^{\circ}$)', '$K$ (MPa)', '$K(T)$', '', 'log');
    else
        fig.ax = [];
    end

    fit1 = draw_exp(fig.ax, T, K, T < T_b, '+', 'exp 1', getMyColor(1));
    fit2 = draw_exp(fig.ax, T, K, T >= T_b, 'o', 'exp 2', getMyColor(2));

    K1 = get_K_interp_local(T1, T_b, fit1, fit2);
end

function logfit = draw_exp(ax, T_all, K_all, ind, marker, tit, clr)
    T = T_all(ind);
    K = K_all(ind);
    
    logfit = polyfit(T, log(K), 1);
    if(~isempty(ax))
        plot(ax, T, K, marker, ...
            'DisplayName', tit, 'Color', clr, 'LineWidth', 1.5);
        plot(ax, T, exp(polyval(logfit, T)), ...
            'Color', clr, 'HandleVisibility', 'off');
    end
end

function K = get_K_interp_local(T, T_b, fit1, fit2)
    if(T < T_b)
        K = exp(polyval(fit1, T));
    else
        K = exp(polyval(fit2, T));
    end
end