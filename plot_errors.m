function plot_errors(Linf,L1,L2,h)
figure
hold on
plot(h,Linf);
plot(h,L1);
plot(h,L2);
xlabel('h');
ylabel('Error');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log');
fontsize(gca,16,"points");
legend("L_{\infty}","L_1","L_2");
pinf=log(Linf(1:end-1)./Linf(2:end))./log(h(1:end-1)./h(2:end));
disp(pinf);
p1=log(L1(1:end-1)./L1(2:end))./log(h(1:end-1)./h(2:end));
disp(p1);
p2=log(L2(1:end-1)./L2(2:end))./log(h(1:end-1)./h(2:end));
disp(p2);