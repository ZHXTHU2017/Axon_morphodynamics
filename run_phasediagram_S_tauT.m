clc;
clear;
close all;

num = 50;
tauvs = logspace(-3, 2, num); 
eta = 2;
Ss = linspace(1 / eta, 7, num); 

real_parts_matrix = nan(num);
imag_parts_matrix = nan(num);

fprintf('Start solving characteristic equation...\n');

options = optimoptions('fsolve', 'Display', 'off');

char_eq = @(lambda_vec, S, tauF, eta) [
    real((lambda_vec(1) + 1i*lambda_vec(2) + 1)^2 + eta * S * exp(-(lambda_vec(1) + 1i*lambda_vec(2)) * (1 + eta) * tauF));
    imag((lambda_vec(1) + 1i*lambda_vec(2) + 1)^2 + eta * S * exp(-(lambda_vec(1) + 1i*lambda_vec(2)) * (1 + eta) * tauF))
];

total_iterations = num^2;
current_iteration = 0;

for i = 1:length(Ss)
    S_val = Ss(i);
    for j = 1:length(tauvs)
        tauF_val = tauvs(j);
        
        current_iteration = current_iteration + 1;
        if mod(current_iteration, 50) == 0
            fprintf('Progress: %d / %d (%.1f%%)\n', current_iteration, total_iterations, (current_iteration/total_iterations)*100);
        end
        
        current_roots = [];
        
        for x_guess = -2:1:2
            for y_guess = -2:1:2
                initial_guess = [x_guess, y_guess];
                [lambda, ~, exitflag] = fsolve(@(lambda) char_eq(lambda, S_val, tauF_val, eta), initial_guess, options);
                
                if exitflag > 0
                    sol = lambda(1) + 1i * lambda(2);
                    if isempty(current_roots) || min(abs(current_roots - sol)) > 1e-4
                        current_roots = [current_roots; sol];
                    end
                end
            end
        end
        
        if ~isempty(current_roots)
            [max_real_part, Index] = max(real(current_roots));
            real_parts_matrix(i, j) = max_real_part;
            imag_parts_matrix(i, j) = imag(current_roots(Index));
        end
    end
end

fprintf('Computation finished.\n');

wc = sqrt(eta * Ss - 1);
tauc = (2./wc/(1+eta)) .* atan(1 ./ wc);

while any(tauc < 0)
    tauc(tauc < 0) = tauc(tauc < 0) + 2*pi./wc(tauc < 0);
end

tauc(tauc == inf) = 1e6;

fprintf('Plotting phase diagram...\n');

stability_map = real_parts_matrix > 0;
display_matrix = abs(imag_parts_matrix / (2*pi));
display_matrix(~stability_map) = NaN;

stable_color = [236,231,242]/255; 

mycolorpoint_rgb = [[33,102,172];
                   [103,169,207];
                   [209,229,240];
                   [253,219,199];
                   [239,138,98];
                   [178,24,43]]; 

mycolorposition = [1 17 26 39 51 64];

mycolormap_r = interp1(mycolorposition, mycolorpoint_rgb(:,1), 1:64, 'linear', 'extrap');
mycolormap_g = interp1(mycolorposition, mycolorpoint_rgb(:,2), 1:64, 'linear', 'extrap');
mycolormap_b = interp1(mycolorposition, mycolorpoint_rgb(:,3), 1:64, 'linear', 'extrap');

mycolor = [mycolormap_r', mycolormap_g', mycolormap_b']/255;
mycolor = round(mycolor*1e4)/1e4;
oscillation_colormap = mycolor;

figure;
ax = gca;
ax.Color = stable_color;
hold on;

imagesc(tauvs, Ss, display_matrix);

colormap(ax, oscillation_colormap);

set(gca, 'XScale', 'log');
axis tight;
box on;

plot(tauc, Ss, 'w-', 'LineWidth', 2);

colorbar;

xlim([min(tauvs) max(tauvs)]);
ylim([0.5 7]);

fprintf('Phase diagram completed.\n');

print('PhaseDiagram_S_tauv_DelayFeedback.eps', '-depsc');


[Tau_grid, S_grid] = meshgrid(tauvs, Ss);
freq_matrix = abs(imag_parts_matrix / (2*pi));

phase_data.S = S_grid;
phase_data.tauT = Tau_grid;
phase_data.frequency = freq_matrix;
phase_data.stability = stability_map;

save('PhaseDiagram_S_tauT_data.mat', '-struct', 'phase_data');

fprintf('Data saved to PhaseDiagram_S_tauT_data.mat\n');
