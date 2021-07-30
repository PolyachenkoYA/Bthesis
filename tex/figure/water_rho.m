clear; close all;

TK2C = 273.15;
T0 = [35, 26, 18, 10, 2];
w_extra_sat = [510, 1000, 1000, 1000, 1000];
K_exp_my = [110, 164, 190, 300, 1290];
do_ind = [0, 1, 1, 1, 1];

eta = 1;
do_rho = 1;
do_L = 0;

%w_extra = [-40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]' + w_extra0;
%w_extra = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610]' + w_extra0;
w_extra_data = {[(0:50:450)'; ([-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]' + (1014 - 180))], ...
                [0, 150, 300, 450, 600, 750]', ...
                [(0:150:450), 650, 800]', ...
                [0, 150, 300, 450, 600, 750]', ...
                [0, 300, 600, 750]'};

%rho = [6.96, 6.36, 5.81, 7.27, 6.00, 7.70, 7.0, 8.0, 7.0, 7.86, 10, 7.3, 7.33, 8.0, 7.66, 8.43]';
%rho = [6.5, 9, 8, 6.1, 7, 10, 8, 7.85, 7, 8.5, 7.9, 8]';
%rho = [5.76, 5.91, 6.77, 6.66, 7.38, 6.65, 6.4, 7.06, 7.17, 7.78]';
rho_data = {[5.76, 5.91, 6.77, 6.66, 7.38, 6.65, 6.4, 7.06, 7.17, 7.78, 6.5, 9, 8, 6.1, 7, 10, 8, 7.85, 7, 8.5, 7.9, 8, 6.96, 6.36, 5.81, 7.27, 6.00, 7.70, 7.0, 8.0, 7.0, 7.86, 10, 7.3, 7.33, 8.0, 7.66, 8.43]', ...
            [3.35, 3.53, 3.58, 3.86, 4.11, 3.86]', ...
            [1.55, 1.89, 2.08, 2.16, 2.58, 2.58]', ... % 2.58
            [0.74, 0.96, 1.15, 1.09, 1.32, 1.18]', ...
            [0.422, 0.48, 0.57, 0.567]'};

%d_rho = [0.456, 0.348, 0.144, 0.800, 1.0, 1.0, 0.1, 1.5, 1.5, 1.18, 3.5, 1.5, 0.942, 0.4, 1.0, 0.182]';
%d_rho = [0.6, 1, 2, 0.3, 1, 3, 2, 0.5, 1, 1.5, 0.6, 2]';
%d_rho = [0.06, 0.06, 0.07, 0.065, 0.07, 0.07, 0.06, 0.07, 0.07, 0.07]';
d_rho_data = {[0.6, 1, 1, 1.2, 0.8, 0.7, 0.56, 1, 2, 1.5, 0.6, 1, 2, 0.3, 1, 3, 2, 0.5, 1, 1.5, 0.6, 2, 0.456, 0.348, 0.144, 0.800, 1.0, 1.0, 0.1, 1.5, 1.5, 1.18, 3.5, 1.5, 0.942, 0.4, 1.0, 0.182]', ...
              [0.08, 0.083, 0.084, 0.087, 0.092, 0.091]', ...
              [0.03, 0.03, 0.03, 0.03, 0.05, 0.05]', ...
              [0.022, 0.025, 0.028, 0.028, 0.031, 0.035]', ...
              [0.018, 0.02, 0.03, 0.03]'};
%d_rho = ones(size(rho)) * 0.4;

if(do_rho)
    fig_rho_Kfit = getFig('$(v_2^{-1/3} - 1) \cdot (5/3 \cdot v_2^{-1/3} - 1) / v_2^2$', '$(\ln(P / P_0) - \ln(1 - v_2) - v_2) / v_2^2$', '$y(x)$');
end
if(do_L)
    fig_Lxy_w = getFig('$w_{extra}$ ($H_2 O$ / cell)', '$L$ (nm)', '$L(w_{extra})$');
    fig_Lxy_rho = getFig('$\rho_{vapor}$ (g/$m^3$)', '$L$ (nm)', '$L(w_{extra})$');
end

N_models = length(T0);
for im = 1:N_models
    if(do_ind(im))
        T = T0(im) + TK2C;
        rho_sat = get_rho_sat(get_opt_rho_sat([]), T) / 1.5;
        %rho_sat = 2.6;    
        not_sat_ind = w_extra_data{im} < w_extra_sat(im);
        w_extra = w_extra_data{im}(not_sat_ind);
        rho = rho_data{im}(not_sat_ind);
        d_rho = d_rho_data{im}(not_sat_ind);

        Ldata = [[7.3304770548, 0.015432846612721736, 7.1722916656, 0.018129885433487338]; ...
                [7.181661062399999, 0.018338095855549075, 7.1532357532, 0.01483626308301014]; ...
                [7.154355831600001, 0.024599868982025933, 7.1407179884, 0.019213137300676995]; ...
                [7.203556148800001, 0.01737954576049267, 7.157062537200001, 0.019548653034002526]; ...
                [7.14633481, 0.016616372016788138, 7.199240510800001, 0.01718012919565751]; ...
                [7.172980352399999, 0.01559303956660839, 7.1008510712, 0.018091184556565953]; ...
                [7.137558040800001, 0.03651943615975383, 7.1998307732, 0.018970557328981178]; ...
                [7.0711723064, 0.018759809651308277, 7.1621753956, 0.01655031349844167]; ...
                [7.072170213600001, 0.017807973590410975, 7.178682769999999, 0.020135832643014787]; ...
                [7.1689189096, 0.016244877394693623, 7.1224542972, 0.016652257547109712]];
        Lx = Ldata(:, 1);
        d_Lx = Ldata(:, 2);
        Ly = Ldata(:, 3);
        d_Ly = Ldata(:, 4);

        d_rho_min = 0.75;

        wrong_ids = d_rho < d_rho_min;
        wght = 1 ./ d_rho;
        Nw = (w_extra + 207 * 8) * 2;
        Na = 6.02 * 1e23;
        Nprot = 8 * 2;
        Vvac = 7.1 * 7.1 * 25 * 1e-27;
        Vprot = 7.1 * 7.1 * (32.5 - 25) * 1e-27;
        h = (Nw - rho * Vvac / 18 * Na) * 18 / (14530 * Nprot);
        d_h = d_rho *  Vvac * Na / (14530 * Nprot);
        %rho_sat = get_opt_rho_sat_2(h, rho, d_rho, 1);
        v2 = 1 ./ (1 + h * 1.38);
        d_v2 = 1.38 * d_h ./ (1 + h * 1.38).^2;
        v1 = 1 - v2;
        %x = Nw / Nprot * (1 / v1 - 1);
        %V1 = Vprot / (Nw + x * Nprot);   % = Vprot / Nw * v1
        V1 = 0.018 / (Na * 1000);   % rho_w = 1 g/cm^3
        phi = rho / rho_sat;
        d_phi = d_rho / rho_sat;
        Kfit_x = (v2 .^ (-1/3) - 1) .* (5/3 * v2 .^ (-1/3) - 1) ./ v2.^2;
        d_Kfit_x = abs((2*(20-28*v2.^(1/3)+9*v2.^(2/3)))./(9*v2.^(11/3))) .* d_v2;
        Kfit_y = (log(phi) - log(1 - v2) * eta - v2 * eta) ./ v2.^2;
        d_Kfit_y = sqrt((((((-2+v2).*v2*eta)./(-1+v2)-2*log(phi./(1-v2)))./v2.^3) .* d_v2) .^ 2 + (d_phi ./ phi ./ v2.^2).^2);
        Kfit_lin = polyfit(Kfit_x, Kfit_y, 1);
        Kfit_K = Kfit_lin(1) / (V1 / (2 * T * 8.31 / Na));
        [~, d_Kfit_lin, Kfit_lin_CI] = ...
            lin_analyze_ax([], Kfit_x, Kfit_y, {'linear', 'linear'}, 1, d_Kfit_y);
        d_Kfit_K = Kfit_K * d_Kfit_lin(1) / Kfit_lin(1);

        fit_obj = fit(w_extra, rho, fittype('poly1'), 'Weight', 1 ./ d_rho);
        linfit = coeffvalues(fit_obj);
        d_linfit = confint(fit_obj);
        d_linfit = (d_linfit(2, :) - d_linfit(1, :)) / 2;
        d_rho_fnc = @(w)sqrt(d_linfit(2) ^ 2 + (d_linfit(1) * w) .^ 2);
        w_ext_fnc = @(rho)([(rho - linfit(2)) / linfit(1),...
                            sqrt((rho * d_linfit(1) / linfit(1)^2).^2 + (d_linfit(2) / linfit(1))^2 + (d_linfit(1) * linfit(2) / linfit(1)^2)^2)]);

        if(do_rho)
            units_K = 1e-6;
            clr = getMyColor(im+2);
            %tit = ['T = ' num2str(T - TK2C) '$C^{\circ}$'];
            %title(fig_rho_Kfit.ax, tit);
                       
            errorbar(fig_rho_Kfit.ax, Kfit_x, Kfit_y, d_Kfit_y, d_Kfit_y, d_Kfit_x, d_Kfit_x, 'o', ...
                'HandleVisibility', 'off', 'Color', clr, 'LineWidth', 1.5);
%             plot(fig_rho_Kfit.ax, Kfit_x, polyval(Kfit_lin, Kfit_x), ...
%                 'DisplayName', ['$K = (' num2str(round(Kfit_K * units_K, -2)) ' \pm ' num2str(round(d_Kfit_K * units_K * 0.85, -2)) ')$ MPa'], ...
%                 'Color', clr, 'LineWidth', 1.5);

            x_bounds = [min(Kfit_x), max(Kfit_x)];
            x_mean = mean(x_bounds);
            y_mean = mean(polyval(Kfit_lin, x_bounds));
            K_exp = get_K_interp(T - TK2C, 0) / units_K;
            K_exp = K_exp_my(im) / units_K;
            plot(fig_rho_Kfit.ax, Kfit_x, (Kfit_x - x_mean) * Kfit_lin(1) * K_exp  / Kfit_K + y_mean, '--', ...
                'DisplayName', ['$T = ' num2str(T - TK2C) ' (C^{\circ}); K_{exp} = ' num2str(round(K_exp * units_K, -1)) '$ MPa'], ...
                'Color', clr, 'LineWidth', 1.5);

    %         fig_rho_w = getFig('$h$ ($m_{w} / m_{prot}$)', '$\rho_w$ (g / $m^3$)', '$\rho_w(w_{extra})$');
    %         errorbar(fig_rho_w.ax, h, rho, d_rho, 'o', 'DisplayName', 'data');
    %         plot(fig_rho_w.ax, h, ones(size(w_extra)) * rho_sat, 'DisplayName', '$\rho_{sat}$');
    %         xlim(fig_rho_w.ax, [0, max(h)]);
    %         ylim(fig_rho_w.ax, [0, max(rho)]);

    %         fig_rho_h = getFig('$\rho_w / \rho_{sat}$', '$h$ ($m_{w} / m_{prot}$)', '$h(\varphi)$');
    %         errorbar(fig_rho_h.ax, phi, h, d_h, d_h, d_phi, d_phi, 'o', 'DisplayName', ['T = ' num2str(T0) '$C^{\circ}$']);
    %         xlim(fig_rho_h.ax, [0, max(phi)]);
    %         ylim(fig_rho_h.ax, [0, max(h)]);    
        end

        if(do_L)
            errorbar(fig_Lxy_w.ax, w_extra, Lx, d_Lx, 'o', 'DisplayName', '$L_x$');
            errorbar(fig_Lxy_w.ax, w_extra, Ly, d_Ly, 'o', 'DisplayName', '$L_y$');

            errorbar(fig_Lxy_rho.ax, rho, Lx, d_Lx, d_Lx, d_rho, d_rho, 'o', 'DisplayName', '$L_x$');
            errorbar(fig_Lxy_rho.ax, rho, Ly, d_Ly, d_Ly, d_rho, d_rho, 'o', 'DisplayName', '$L_y$');
        end

        % plot(fig.ax, w_extra, polyval(linfit, w_extra), ...
        %     'DisplayName', 'linfit', 'Color', getMyColor(2));
        % plot(fig.ax, w_extra, polyval(linfit, w_extra) + d_rho_fnc(w_extra), '--', ...
        %     'HandleVisibility', 'off', 'Color', getMyColor(2));
        % plot(fig.ax, w_extra, polyval(linfit, w_extra) - d_rho_fnc(w_extra), '--', ...
        %     'HandleVisibility', 'off', 'Color', getMyColor(2));

        % yt = w_ext_fnc(polyval(linfit, w_extra));
        % plot(fig.ax, yt(:, 1), polyval(linfit, w_extra), ...
        %     'DisplayName', 'test', 'Color', getMyColor(3));
        % plot(fig.ax, yt(:, 1) - yt(:, 2), polyval(linfit, w_extra), ...
        %     'HandleVisibility', 'off', 'Color', getMyColor(3));
        % plot(fig.ax, yt(:, 1) + yt(:, 2), polyval(linfit, w_extra), ...
        %     'HandleVisibility', 'off', 'Color', getMyColor(3));
    end
end

