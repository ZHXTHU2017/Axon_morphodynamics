clc;
clear;
close all;
% 
tau = 100; 

K = 0.5; 
KR = 0.12;
n = 4; % 
ka = 1 / tau; 
kr = 1 / tau; 
vm = 1; 
Lk = 27; 
zeta_r = 5.2; 
Tr = 1; 
zeta_a = 10; 
vg = 0.21; 
ms = 10; 
v0 = 0.1;
eta = 2.0; 
%% Normalization
ka_nondim  = ka * tau;
kr_nondim  = kr * tau; 
v_nondim = vm * tau / Lk;
JR_nondim = 1.0; 
zeta_r_nondim = zeta_r * Lk ^ 2 / (tau * Tr);
zeta_a_nondim = zeta_a * Lk / (tau * Tr);
ms_nondim = ms * K / Tr;
v0_nondim = v0 * tau / Lk;
vg_nondim = vg * K * tau / Lk;
Ja_nondim = 1.6;
Jr_nondim = 2.4;
% Steady state
system_eqs = @(x,Ja,Jr,n,eta) [
    x(1) - eta * Ja * (1 - (x(2)^n) / (1 + x(2)^n));
    x(2) - Jr * (x(1)^n) / (1 + x(1)^n)
];
x0 = [1; 1];
options = optimoptions('fsolve', 'Display', 'off');
fsolve_handle = @(x) system_eqs(x, Ja_nondim, Jr_nondim, n, eta);
[solution, fval, exitflag] = fsolve(fsolve_handle, x0, options);

if exitflag > 0
    ca_sol = solution(1);
    cr_sol = solution(2);
    
    fprintf('------------------------------------\n');
    fprintf('ca* = %.6f\n', ca_sol);
    fprintf('cr* = %.6f\n', cr_sol);
    fprintf('------------------------------------\n');
else
    fprintf('fsolve exitflag: %d\n', exitflag);
end
df = @(c) n * c^(n-1) / (1 + c ^n) ^2;
Pa_nondim = Ja_nondim * df(cr_sol);
Pr_nondim = Jr_nondim * df(ca_sol);
S_nondim = Pa_nondim * Pr_nondim;
fprintf('Dimensionless S = %.4f \n', S_nondim);
fprintf('\n');

params.Ja = Ja_nondim;
params.Jr = Jr_nondim;
params.n = n;       
params.v = v_nondim;      
params.zeta_a = zeta_a_nondim; 
params.zeta_r = zeta_r_nondim; 
params.JR = JR_nondim;     
params.g = ms_nondim;     
params.vg = vg_nondim;     
params.v0 = v0_nondim;     
params.KR = KR/K;
params.eta = eta;

cR_sol = params.JR * (1 - ca_sol^n/(params.KR^n+ca_sol^n));

F = zeta_a_nondim*vg_nondim*ca_sol - ms_nondim*cR_sol*(1-vg_nondim*ca_sol/v0_nondim)-1;

t_start = 0;
t_end = 5000;

l0=0.1;
initial_conditions = [0; 0; 0; l0]; 

history = @(t) initial_conditions;

options = ddeset('InitialY', initial_conditions);
sol = ddesd(@(t,y,Z) ddefun(t,y,Z,params), @(t,y) delays(t,y,params), history, [t_start, t_end], options);

t = sol.x';
ca = sol.y(1,:)';
cr = sol.y(2,:)';
cR = sol.y(3,:)';
L = sol.y(4,:)';
Gamma = params.zeta_a + params.zeta_r ./ tanh(L) + params.g .* cR ./ params.v0;
L_dot = (1./Gamma) .* (params.zeta_a * params.vg .* ca - ...
    params.g .* cR .* (1 - params.vg .* ca / params.v0) - 1);
Vr_nondim = params.vg .* ca - L_dot;
Vp_nondim = params.vg .* ca;
Fa_nondim = params.zeta_a * (params.vg .* ca - L_dot);
T_nondim = params.zeta_r * L_dot .* tanh(L) + 1;
Tm_nondim = params.g .* cR .* (1 - (params.vg .* ca - L_dot) ./ params.v0);
F_nondim = params.zeta_a * params.vg .* ca - ...
    params.g .* cR .* (1 - params.vg .* ca / params.v0) - 1;
% 绘制结果

L_smooth = smooth(t, L, 0.07, 'lowess');

figure
subplot(2,1,1)
plot(sol.x, ca, 'b', 'LineWidth', 2)
hold on;
plot(sol.x, cr, 'r', 'LineWidth', 2)
plot(sol.x, cR, 'g', 'LineWidth', 2)
subplot(2,1,2)
plot(sol.x, L, 'm', 'LineWidth', 2);hold on;
grid off;
% plot(t, L_smooth, 'b', 'LineWidth', 2);

% figure;
% subplot(3,1,1)
% plot(sol.x, Fa_nondim, 'b', 'LineWidth', 2)
% subplot(3,1,2)
% plot(sol.x, T_nondim, 'm', 'LineWidth', 2)
% subplot(3,1,3)
% plot(sol.x, Tm_nondim, 'm', 'LineWidth', 2)
% figure;
% subplot(3,1,1)
% plot(sol.x, L_dot, 'b', 'LineWidth', 2)
% subplot(3,1,2)
% plot(sol.x, Vp_nondim, 'm', 'LineWidth', 2)
% subplot(3,1,3)
% plot(sol.x, Vr_nondim, 'm', 'LineWidth', 2)

dt_uniform = 0.5; 
t_u = (t_start+10:dt_uniform:t_end)';
interp_func = @(data) interp1(t, data, t_u, 'pchip'); 

cr_u = interp_func(cr);
ca_u = interp_func(ca);
cR_u = interp_func(cR);
L_u  = interp_func(L_smooth);
L_dot_u = interp_func(L_dot);
Vr_u = interp_func(Vr_nondim);
Vp_u = interp_func(Vp_nondim);
Fa_u = interp_func(Fa_nondim);
T_u  = interp_func(T_nondim);
F_u  = interp_func(F_nondim);
Tm_u = interp_func(Tm_nondim);
Gamma_u = interp_func(Gamma);

avg_ca = get_trend_exact(t_u, ca_u);
avg_cr = get_trend_exact(t_u, cr_u);
avg_cR = get_trend_exact(t_u, cR_u);
avg_L = L_u;
avg_L_dot = get_trend_exact(t_u, L_dot_u);
avg_Gamma = get_trend_exact(t_u, Gamma_u);

avg_Vr = get_trend_exact(t_u, Vr_u);
avg_Vp = get_trend_exact(t_u, Vp_u);
avg_Fa = get_trend_exact(t_u, Fa_u);
avg_T  = get_trend_exact(t_u, T_u); 
avg_F  = get_trend_exact(t_u, F_u);
avg_Tm = get_trend_exact(t_u, Tm_u);

avg_Fa_smooth = smooth(t_u, avg_Fa, 0.2, 'lowess');
avg_Tm_smooth = smooth(t_u, avg_Tm, 0.2, 'lowess');


R_dFa = gradient(avg_Fa_smooth) ./ gradient(t_u);
R_dTm = gradient(avg_Tm_smooth) ./ gradient(t_u);



global_ptp = max(cr_u) - min(cr_u);
prom_threshold = 0.05 * global_ptp; 
if prom_threshold < 1e-6, prom_threshold = 1e-6; end

[~, locs_cr] = findpeaks(cr_u, 'MinPeakProminence', prom_threshold);


time_freq_seq = t_u;
freq_cr_seq = zeros(size(t_u));

if length(locs_cr) >= 2
    
    peak_times = t_u(locs_cr);                 
    periods = diff(peak_times);               
    inst_freqs = 1 ./ periods;             
    
    mid_times = peak_times(1:end-1) + periods / 2;
    
    freq_cr_seq(t_u < peak_times(1)) = 0;
    
    if length(mid_times) > 1
        
        freq_cr_seq_interp = interp1(mid_times, inst_freqs, t_u, 'pchip');
        osc_idx = (t_u >= peak_times(1)) & (t_u <= peak_times(end));
        freq_cr_seq(osc_idx) = freq_cr_seq_interp(osc_idx);   
        freq_cr_seq(t_u > peak_times(end)) = inst_freqs(end);
    else
        
        freq_cr_seq(t_u >= peak_times(1)) = inst_freqs(1);
    end
else
    fprintf('警告：未检测到足够的有效振荡周期，频率保持为 0。\n');
end

freq_cr_seq = smoothdata(freq_cr_seq, 'gaussian', 30);


%% 3. Plot
plot(avg_L, avg_T, 'LineWidth', 2.5, 'Color', [0.6 0 0.8]); 
xlabel('Average Length (L)'); 
ylabel('Average Tension (T)');
grid off;box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% figure
% subplot(2,1,1);
% plot(t_u, cr_u, 'Color', [0.8 0.8 0.8]); hold on; % 原始
% plot(t_u, avg_cr, 'b-', 'LineWidth', 2);          % 平均
% ylabel('c_r'); title('浓度 cr');
% grid on; 
% 
% subplot(2,1,2);
% plot(t_u, T_u, 'Color', [0.8 0.8 0.8], 'LineWidth', 1); hold on;  
% plot(t_u, avg_T, 'r-', 'LineWidth', 2);   
% ylabel('Tension');
% xlabel('Time'); 
% xlim([0, max(t_u)]); 
% grid on;

figure;
% subplot(2,1,1);
% plot(t_u, avg_L, 'b-', 'LineWidth', 2);  
% title("平均长度")
% grid on; 
% subplot(2,1,2);
yyaxis left; plot(t_u, avg_T, 'r-', 'LineWidth', 2);   ylabel('Mean T');
yyaxis right; plot(t_u, avg_L_dot, 'r-', 'LineWidth', 2);  ylabel('Mean dL/dt');
xlabel('Time'); grid off;


% figure;
% subplot(2,1,1);
% plot(t_u, ca_u, 'Color', [0.8 0.8 0.8], 'LineWidth', 2)
% hold on;
% plot(t_u, avg_ca, 'r', 'LineWidth', 2)
% subplot(2,1,2);
% plot(t_u, cR_u, 'Color', [0.8 0.8 0.8], 'LineWidth', 2)
% hold on;
% plot(t_u, avg_cR, 'r', 'LineWidth', 2)

% figure;
% yyaxis left; plot(avg_L, avg_Fa, '-b', 'LineWidth', 1); ylabel('Mean Value Fa');
% hold on; plot(avg_L, avg_T, 'k-', 'LineWidth', 1);
% yyaxis right; plot(avg_L, avg_Tm, '-r', 'LineWidth', 1); ylabel('Mean Value Tm');
% legend('Mean Fa', 'Mean T', 'Mean Tm');grid on;

figure;
yyaxis left; scatter(freq_cr_seq, avg_ca, 'b'); ylabel('Mean ca');
yyaxis right; scatter(freq_cr_seq, avg_cR, 'r'); ylabel('Mean cR');
legend('Mean ca', 'Mean cR');grid off;

% figure;
% subplot(3,1,1)
% plot(t_u, L_dot_u, 'Color', [0.8 0.8 0.8], 'LineWidth', 2)
% hold on;
% plot(t_u, avg_L_dot, 'r', 'LineWidth', 2)
% 
% subplot(3,1,2)
% plot(t_u, Vr_u,'Color', [0.8 0.8 0.8], 'LineWidth', 2)
% hold on;
% plot(t_u, avg_Vr, 'r', 'LineWidth', 2)
% 
% subplot(3,1,3)
% plot(t_u, Vp_u, 'Color', [0.8 0.8 0.8], 'LineWidth', 2)
% hold on;
% plot(t_u, avg_Vp, 'r', 'LineWidth', 2)


figure;
plot(t_u, R_dFa,'r', 'LineWidth', 2);hold on;
plot(t_u, R_dTm,'b', 'LineWidth', 2);hold off;
grid off;
figure
% subplot(2,1,1);
% plot(t_u, cr_u, 'Color', [0.7 0.7 0.7]); hold on;
% ylabel('cr');
% xlim([t_start, t_end]); grid on;
% 
% subplot(2,1,2);
yyaxis left; plot(time_freq_seq, freq_cr_seq, 'b-', 'LineWidth', 2);
% start_osc_idx = find(freq_cr_seq > 0.01, 1);
% if ~isempty(start_osc_idx)
%     xline(time_freq_seq(start_osc_idx), 'r--', 'Label', 'Oscillation initiation');
% end
ylabel('Frequency (Hz)'); xlabel('Time');
ylim([0, max(freq_cr_seq)*1.2]);

yyaxis right; plot(t_u, avg_L, 'r-', 'LineWidth', 2);ylabel('Mean length, L');
xlim([t_start, t_end]); 
grid off;

function trend = get_trend_exact(t, y)

    prom = 0.001 * (max(y) - min(y)); 
    if prom < 1e-6, prom = 1e-6; end
    [~, locs] = findpeaks(y, 'MinPeakProminence', prom);
    

    if length(locs) < 5
        trend = smoothdata(y, 'gaussian', 50);
        return;
    end


    num_cycles = length(locs)-1;
    cycle_time_points = zeros(num_cycles, 1);
    cycle_raw_means = zeros(num_cycles, 1);
    
    for i = 1:num_cycles
        idx1 = locs(i);
        idx2 = locs(i+1);
        t_seg = t(idx1:idx2);
        y_seg = y(idx1:idx2);
        
        dt_cycle = t_seg(end) - t_seg(1);
        area = trapz(t_seg, y_seg);
        
        cycle_raw_means(i) = area / dt_cycle;
        cycle_time_points(i) = mean(t_seg); 
    end
    
    window_span = 15; 
    if length(cycle_raw_means) > window_span
        cycle_smooth_means = smoothdata(cycle_raw_means, 'gaussian', window_span);
    else
        cycle_smooth_means = cycle_raw_means;
    end
    
    first_cycle_idx = locs(1);
    t_split = t(first_cycle_idx);

    t_head = t(t <= t_split);
    y_head = y(t <= t_split);
    
   
    y_head = smoothdata(y_head, 'gaussian', 5);
    

    t_source = [t_head; cycle_time_points];
    y_source = [y_head; cycle_smooth_means];
    
    [t_source, sort_idx] = sort(t_source);
    y_source = y_source(sort_idx);
    [t_source, unique_idx] = unique(t_source);
    y_source = y_source(unique_idx);

    if t_source(end) < t(end)
        t_source = [t_source; t(end)];
        y_source = [y_source; y_source(end)];
    end

    try

        trend = interp1(t_source, y_source, t, 'pchip');
    catch
        trend = smoothdata(y, 'gaussian', 50);
    end
    
    trend = smoothdata(trend, 'gaussian', 20); 
end

function dydt = ddefun(t, y, Z, params)

ca = y(1);
cr = y(2);
cR = y(3);
L = y(4);
cr_delayed =  Z(2,1);
ca_delayed =  Z(1,2);

dca_dt = -ca + params.eta * params.Ja / (1 + cr_delayed^params.n);

dcr_dt = -cr + params.Jr * (ca_delayed^params.n) / (1 + ca_delayed^params.n);
dcR_dt = -cR + params.JR * (1 - ca^params.n / (params.KR^params.n + ca^params.n));

Gamma = params.zeta_a + params.zeta_r / tanh(L) + params.g * cR / params.v0;

dL_dt = (1/Gamma) * (params.zeta_a * params.vg * ca - ...
    params.g * cR * (1 - params.vg * ca / params.v0) - 1);
dydt = [dca_dt; dcr_dt; dcR_dt; dL_dt];
end
function d = delays(t, y, params)
L = y(4);
d = [t-params.eta*L / params.v
    t-L / params.v];
end