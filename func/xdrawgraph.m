function xdrawgraph(xs,yrange,method_list,field,ti,lx,ly)

box('on');
hold('all');

p= zeros(size(method_list));
for i= 1:length(method_list)
    p(i)= plot(xs,method_list(i).(field),method_list(i).marker,...
        'displayname',method_list(i).name);
end
ylim(yrange);
xlim(xs([1 end]));
set(gca,'xtick',xs)

title(ti,'FontSize',12,'FontName','Arial');
xlabel(lx);
ylabel(ly);
legend(p);

return
