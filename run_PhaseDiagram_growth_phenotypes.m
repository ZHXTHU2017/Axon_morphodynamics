clc;
clear;
close all;

tau = 100;    
K = 0.5;       
KR = 0.12;
n = 4;         
vm = 1;        
Lk = 27;       
zeta_r = 5.2;  % frictional viscosity, kPa s = nN um-2 s
Tr = 1;        % rest tension, nN
zeta_a = 10;   % Adhesive frictional viscosity，pN s nm-1
vg = 0.21;     % Signal-modulated velocity constant，um s-1 / (n.u um-1)
ms = 10;       % Signal-modulated contractility constant，nN / (n.u um-1)
v0 = 0.1;      % Unloaded retrograde flow velocity , um s-1
eta = 2.0;    


v_nondim = vm * tau / Lk;
JR_nondim = 1.0; 
Ja_nondim = 1.6; 
zeta_r_nondim = zeta_r * Lk ^ 2 / (tau * Tr);
zeta_a_nondim = zeta_a * Lk / (tau * Tr);
ms_nondim = ms * K / Tr;
v0_nondim = v0 * tau / Lk;
vg_nondim = vg * K * tau / Lk;

S_threshold = 1 / eta; 

N_points = 50; 
Jr_list = logspace(-1, 2, N_points);       
tauv_list = logspace(-3, 2, N_points);     

fhill = @(c) c.^n ./ (1 + c.^n);
df = @(c) n * c.^(n-1) ./ (1 + c.^n).^2;
fsolve_opts = optimoptions('fsolve', 'Display', 'off');

S_list = zeros(1, N_points);
Fstable_list = zeros(1, N_points);

for j = 1:N_points
    Jr_curr = Jr_list(j);
    system_eqs = @(x) [
        x(1) - eta * Ja_nondim * (1 - fhill(x(2)));
        x(2) - Jr_curr * fhill(x(1))
    ];
    [sol_ss, ~, exitflag] = fsolve(system_eqs, [1; 1], fsolve_opts);
    if exitflag > 0
        ca_ss = sol_ss(1);
        cr_ss = sol_ss(2);

        Pa_nondim = Ja_nondim * df(cr_ss);
        Pr_nondim = Jr_curr * df(ca_ss);
        S_list(j) = Pa_nondim * Pr_nondim;
        
        cR_ss = JR_nondim * (1 - ca_ss^n / ((KR/K)^n + ca_ss^n));
        Fstable_list(j) = zeta_a_nondim * vg_nondim * ca_ss - ms_nondim * cR_ss * (1 - vg_nondim * ca_ss / v0_nondim) - 1;
    else
        S_list(j) = NaN;
        Fstable_list(j) = NaN;
    end
end

max_S = max(S_list); 

[tauv_grid, Jr_grid] = meshgrid(tauv_list, Jr_list);
[~, S_grid] = meshgrid(tauv_list, S_list); 
[~, Fstable_grid] = meshgrid(tauv_list, Fstable_list); 

total_iter = N_points * N_points;

mean_ca_vec   = NaN(total_iter, 1);
mean_cR_vec   = NaN(total_iter, 1);
mean_Fnet_vec = NaN(total_iter, 1);

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;
end
D = parallel.pool.DataQueue;
h_wait = waitbar(0, '');
afterEach(D, @(~) update_waitbar(h_wait, total_iter));

parfor k = 1:total_iter
    Jr_curr = Jr_grid(k);
    tauv_curr = tauv_grid(k);
    S_curr = S_grid(k); 
    
    lags = tauv_curr;
    lag_vec = [lags, eta * lags];
    history = [0; 0; 0];

    ca_cr_func = @(t,c,Z) [ ...
        -c(1) + eta * Ja_nondim*(1 - Z(2,2)^n / (1 + Z(2,2)^n)); ...
        -c(2) + Jr_curr*(Z(1,1)^n / (1 + Z(1,1)^n));
        -c(3) + JR_nondim * (1 - c(1)^n / ((KR/K)^n + c(1)^n))];

    t_end = max(500, 200 * tauv_curr);
    opts = ddeset('RelTol', 1e-4, 'AbsTol', 1e-4);

    sol = dde23(ca_cr_func, lag_vec, history, [0, t_end], opts);

  
    t_start_ana = t_end * 0.7;
    num_samples = 2000;
    t_uniform = linspace(t_start_ana, t_end, num_samples);
    sol_uniform = deval(sol, t_uniform);

    ca_data = sol_uniform(1, :);
    cR_data = sol_uniform(3, :);

    Fnet_data = zeta_a_nondim * vg_nondim * ca_data - ...
        ms_nondim * cR_data .* (1 - vg_nondim * ca_data / v0_nondim) - 1;


    mean_ca_vec(k)   = mean(ca_data);
    mean_cR_vec(k)   = mean(cR_data);
    mean_Fnet_vec(k) = mean(Fnet_data);

    send(D, 1);
end

if isvalid(h_wait)
    close(h_wait);
end
fprintf('迭代计算完成！正在生成云图...\n');

num = 50;
Ss= linspace(0.5, 8, num); 


wc = sqrt(eta * Ss - 1);
tauc = (2./wc/(1+eta)) .* atan(1 ./ wc); 

while tauc < 0
    tauc = tauc + 2*pi/wc;
end

tauc(tauc == inf)= 1e6;

% mean_ca_mat   = reshape(mean_ca_vec, [N_points, N_points]);
% mean_cR_mat   = reshape(mean_cR_vec, [N_points, N_points]);
mean_Fnet_mat = reshape(mean_Fnet_vec, [N_points, N_points]);

% figure('Name', '2D Parameter Sweep (tau_v vs S)', 'Position', [100, 100, 1500, 450], 'Color', 'w');
% 
% subplot(1, 2, 1);
% pcolor(tauv_grid, S_grid, mean_ca_mat); 
% shading interp;
% set(gca, 'XScale', 'log'); 
% ylim([S_threshold, max_S]);
% colormap(gca, 'jet');
% cb1 = colorbar; title(cb1, 'Mean c_a');
% xlabel('Time Delay \tau_v (Log Scale)'); 
% ylabel('Feedback Strength \itS');
% 
% set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 11);
% box on;
% hold on;
% 
% plot(tauc, Ss, 'w-','LineWidth', 2)
% 
% subplot(1, 2, 2);
% pcolor(tauv_grid, S_grid, mean_cR_mat); 
% shading interp;
% set(gca, 'XScale', 'log'); 
% ylim([S_threshold, max_S]); 
% colormap(gca, 'jet');
% cb2 = colorbar; title(cb2, 'Mean c_R');
% xlabel('Time Delay \tau_v (Log Scale)'); 
% ylabel('Feedback Strength \itS');
% 
% set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 11);
% box on;
% hold on;
% 
% plot(tauc, Ss, 'w-','LineWidth', 2)


figure;
pcolor(tauv_grid, S_grid, mean_Fnet_mat); 
shading interp;
set(gca, 'XScale', 'log'); 
ylim([S_threshold, max_S]);

max_F = max(abs(mean_Fnet_mat(:)), [], 'omitnan');
if isnan(max_F) || max_F == 0, max_F = 1; end
clim([-max_F, max_F]); 

cmap_div = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); 
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)']; 
colormap(gca, cmap_div);

cb3 = colorbar; title(cb3, 'Net Force');
xlabel('Time Delay \tau_v (Log Scale)'); 
ylabel('Feedback Strength \itS');
set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 11);
box on;

hold on;
num_smooth = 300;
tauv_smooth = logspace(log10(min(tauv_list)), log10(max(tauv_list)), num_smooth);
S_smooth    = linspace(min(S_list), max(S_list), num_smooth);
[tauv_grid_sm, S_grid_sm] = meshgrid(tauv_smooth, S_smooth);

Fnet_interp = interp2(tauv_grid, S_grid, mean_Fnet_mat, tauv_grid_sm, S_grid_sm, 'spline');

Fnet_smoothed = imgaussfilt(Fnet_interp, 2); 

[~, h_dyn] = contour(tauv_grid_sm, S_grid_sm, Fnet_smoothed, [0 0], 'k--', 'LineWidth', 2); 

[~, h_ss] = contour(tauv_grid, S_grid, Fstable_grid, [0 0], 'r-.', 'LineWidth', 2); 


plot(tauc, Ss, 'b-','LineWidth', 2)


L0_list = tauv_list * v_nondim;


valid_idx = ~isnan(S_list);
S_valid = S_list(valid_idx);
Fnet_valid = mean_Fnet_mat(valid_idx, :);
Fstable_valid = Fstable_grid(valid_idx, :);

[S_clean, unique_idx] = unique(S_valid);
Fnet_clean = Fnet_valid(unique_idx, :);
Fstable_clean = Fstable_valid(unique_idx, :);


num_fine = 300; 
S_fine = linspace(min(S_clean), max(S_clean), num_fine);
L0_fine = logspace(log10(min(L0_list)), log10(max(L0_list)), num_fine);
[S_grid_fine, L0_grid_fine] = meshgrid(S_fine, L0_fine);

Fnet_fine = interp2(L0_list, S_clean, Fnet_clean, L0_grid_fine, S_grid_fine, 'spline');
Fstable_fine = interp2(L0_list, S_clean, Fstable_clean, L0_grid_fine, S_grid_fine, 'spline');

C_fnet = contourc(S_fine, L0_fine, Fnet_fine, [0 0]);

idx = 1;
left_S_pts = []; left_L0_pts = [];
right_S_pts = []; right_L0_pts = [];
while idx < size(C_fnet, 2)
    num_pts = C_fnet(2, idx);
    S_vals = C_fnet(1, idx+1 : idx+num_pts);
    L0_vals = C_fnet(2, idx+1 : idx+num_pts);
  
    if mean(S_vals) < 2.5 
        left_S_pts = [left_S_pts, S_vals];
        left_L0_pts = [left_L0_pts, L0_vals];
    else                   
        right_S_pts = [right_S_pts, S_vals];
        right_L0_pts = [right_L0_pts, L0_vals];
    end
    idx = idx + num_pts + 1;
end

if ~isempty(left_S_pts)
    Scr2 = min(left_S_pts); 
else
    Scr2 = 1.9; 
end

[~, min_idx] = min(abs(Fstable_list));
Scr1 = S_list(min_idx);


if ~isempty(right_S_pts)
    [Scr3, max_idx] = max(right_S_pts);
    L0_peak = right_L0_pts(max_idx); 
   
    valid_curve_idx = right_S_pts >= Scr1;
    if any(valid_curve_idx)
        L0_int = min(right_L0_pts(valid_curve_idx));
    else
        L0_int = min(L0_fine);
    end

    [right_L0_unq, unq_idx] = unique(right_L0_pts);
    right_S_unq = right_S_pts(unq_idx);
else
    
    Scr3 = 4.1;
    L0_peak = max(L0_fine);
    L0_int = min(L0_fine);
    right_L0_unq = [L0_int, L0_peak];
    right_S_unq = [Scr1, Scr3];
end


phase_map_fine = zeros(size(Fnet_fine));

for i = 1:size(S_grid_fine, 1)      
    for j = 1:size(S_grid_fine, 2)  
        s_val = S_grid_fine(i, j);
        l0_val = L0_grid_fine(i, j);
        
        
        if l0_val <= L0_int
            S_bound = Scr1; % 
        elseif l0_val >= L0_peak
            S_bound = Scr3; 
        else
            
            s_curve_val = interp1(right_L0_unq, right_S_unq, l0_val, 'linear', NaN);
            if isnan(s_curve_val) || s_curve_val < Scr1
                S_bound = Scr1;
            else
                S_bound = s_curve_val;
            end
        end

        if s_val < Scr2
            phase_map_fine(i, j) = 3; 
        elseif s_val > S_bound
            phase_map_fine(i, j) = 1; 
        else
            phase_map_fine(i, j) = 2; 
        end
    end
end

figure('Name', 'Phase Diagram: Strict 3-Segment Boundary', 'Color', 'w', 'Position', [150, 150, 750, 650]);

custom_cmap = [
    186/255, 143/255, 106/255; % 1: Collapse
    139/255, 106/255, 186/255; % 2: Homeostasis 
    102/255, 178/255, 102/255  % 3: Sustained growth
];

% 1. 
pcolor(S_grid_fine, L0_grid_fine, phase_map_fine);
shading flat; 
hold on;
y_lims = [min(L0_fine), max(L0_fine)];

% 2.  
plot([Scr2, Scr2], y_lims, 'r-', 'LineWidth', 4);

plot([Scr1, Scr1], [min(L0_fine), L0_int], 'c-', 'LineWidth', 4);

curve_idx = right_L0_pts >= L0_int & right_L0_pts <= L0_peak;

[L0_plot_sort, sort_i] = sort(right_L0_pts(curve_idx));
S_plot_sort = right_S_pts(find(curve_idx));
S_plot_sort = S_plot_sort(sort_i);
plot(S_plot_sort, L0_plot_sort, 'r-', 'LineWidth', 4);

plot([Scr3, Scr3], [L0_peak, max(L0_fine)], 'r-', 'LineWidth', 4);


colormap(custom_cmap);
set(gca, 'YScale', 'log'); 
xlim([min(S_fine), max(S_fine)]);
ylim(y_lims);

xlabel('Feedback strength, \itS', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Initial length, \itL_0', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1.5, 'FontSize', 13, 'FontWeight', 'bold');
box on;


y_text = median(L0_fine) * 1.5; 
text(Scr2 - 0.4, y_text, 'Sustained growth', 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
text(Scr1 + 0.4, y_text, 'Homeostasis', 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
text(Scr3 + 1.2, y_text, 'Collapse', 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');

y_top = max(L0_fine) * 1.2;
text(Scr2, y_top, '\itS_{cr2}', 'Color', 'r', 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(Scr1, y_top, '\itS_{cr1}', 'Color', 'c', 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(Scr3, y_top, '\itS_{cr3}', 'Color', 'r', 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

set(gca, 'Clipping', 'off');


function update_waitbar(h, total)
    persistent p
    if isempty(p)
        p = 1;
    else
        p = p + 1;
    end
    if isvalid(h)
        waitbar(p/total, h, sprintf('2D %d / %d (%.1f%%)', p, total, p/total*100));
    end
    if p == total
        p = []; 
    end
end
