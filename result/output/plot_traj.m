%
% TUM VI dataset
gt = load('./maj4/ground_truth.tum');
vins_mono = load('./maj4/VINS-Mono.tum');
vins_mono_sd = load('./maj4/VINS-Mono-SD.tum');

plot_xyz(gt(:,2), gt(:, 3), gt(:, 4), 'k'); 
hold on; 
plot_xyz(vins_mono(:,2), vins_mono(:, 3), vins_mono(:, 4), 'g');
hold on; 
plot_xyz(vins_mono_sd(:,2), vins_mono_sd(:, 3), vins_mono_sd(:, 4), 'r');
% grid on; 
legend('groundtruth', 'VINS-Mono', 'VINS-Mono-SD');
title('Majistrale 4');




