

f = load('mean_err_gold_trans_samp.log'); 
s = f(:,1);
re = f(end:-1:1,2);
te = f(end:-1:1,3);
se = f(end:-1:1,4); 

plot(s(1:13), re(1:13), 'k:');
hold on;
plot(s(1:13), te(1:13), 'r-');
hold on;
plot(s(1:13), se(1:13), 'g--');
xlabel('sigma (pixels)');
ylabel('error (pixels)'); 
legend('reprojection error', 'transfer error', 'Sampson error');