

f = load('mean_err_gold_trans_samp.log'); 
s = f(:,1);
re = f(end:-1:1,2);
te = f(end:-1:1,3);
se = f(end:-1:1,4); 

plot(s(1:12), re(1:12), 'k:');
hold on;
plot(s(1:12), te(1:12), 'r-');
hold on;
plot(s(1:12), se(1:12), 'g--');
xlabel('sigma (pixels)');
ylabel('error (pixels)'); 
legend('Reprojection Error', 'Transfer distance', 'Sampson distance');