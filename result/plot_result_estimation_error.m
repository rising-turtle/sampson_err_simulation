
%%
% plot the estimation square error - square distance of the estimated from the
% true value 
% f = load('mean_err_gold_trans_samp.log'); 
% f = load('geometric_dis_mean_err_gold_trans_samp.log'); 
% f = load('mean_error_comparison.log');
f = load('mean_estimation_error_comparison.log');
s = f(:,1);
re = f(:,2);
te = f(:,3);
se = f(:,4); 
% se_g = f(end:-1:1,5); 
se_g = f(:,5); 

figure; 
plot(s(1:12), re(1:12), 'k:');
hold on;
plot(s(1:12), te(1:12), 'r-');
hold on;
plot(s(1:12), se(1:12), 'g--');
hold on; 
%plot(s(1:12), se_g(1:12), 'm-.');
xlabel('sigma (pixels)');
ylabel('mean estimation square error (pixels)'); 
legend('Reprojection Error', 'Transfer Distance', 'Sampson Distance');
%legend('Reprojection Error', 'Transfer Distance', 'Sampson Distance (perspective projection)', ...
%            'Sampson Distance (epipolar)');