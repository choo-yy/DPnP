clear; clc;
addpath epnp;
addpath lhm;
addpath func;
addpath dls_pnp_matlab;
addpath RDLT;
addpath OPnP;
addpath DPnP;

% experimental parameters
nls= 0.5:0.5:4;
npts= 10; npt = npts;
num= 1000;

% compared methods
A= zeros(size(npts));
B= zeros(num,1);
name=   {       'EPnP', 'RDLT+GN', 'DLS', 'LHM', 'RPnP',     'OPnP',        'DPnP',       'DPnP+GN'};
f=      { @EPnP_planar,  @RDLT_GN,  @DLS,  @LHM,  @RPnP,  @OPnP1res,  @DPnP_planar, @DPnP_GN_planar};
marker= {         '-o',      '-x',  '->',  '-+',   '-^',       '-s',          '-d',            '-p'};
color = {"#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#0072BD", "#D95319"};



method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker,'color', color);

% experiments
for i= 1:length(nls)
    
    nl= nls(i);
    fprintf('nl = %d: ',nl);
    
    for j= 1:num
        
        % camera's parameters
        width= 640;
        height= 480;
        f= 800;
        
        % generate 3d coordinates in camera space
        XXw= [xrand(2,npt,[-2 2]); zeros(1,npt)];
		R= rodrigues(randn(3,1));
        t= [rand-0.5;rand-0.5;rand*8+4];
        Xc= R*XXw+repmat(t,1,npt);
        
        % projection
        xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;
		xxn= xxn/f;

        % pose estimation
        for k= 1:length(method_list)
            [R1,t1]= method_list(k).f(XXw,xxn);
            y= cal_pose_err([R1 t1],[R t]);
            method_list(k).r(j)= y(1);
            method_list(k).t(j)= y(2);
        end

        showpercent(j,num);
    end
    fprintf('\n');
    
    % save result
    for k= 1:length(method_list)
        method_list(k).mean_r(i)= mean(method_list(k).r);
        method_list(k).mean_t(i)= mean(method_list(k).t);
        method_list(k).med_r(i)= median(method_list(k).r);
        method_list(k).med_t(i)= median(method_list(k).t);
        method_list(k).std_r(i)= std(method_list(k).r);
        method_list(k).std_t(i)= std(method_list(k).t);
    end
end

%%
close all;
yrange= [0 45];

i= 0; w= 1000; h= 350;
set(groot,'defaultLineLineWidth',2)
figure('color','w','position',[w*i,100,w,h]);i=i+1;
subplot(121)
xdrawgraph(nls,yrange,method_list,'mean_r','Mean Rotation Error',...
    'Number of Points','Rotation Error (degrees)');
subplot(122)
yrange= [0 5];
xdrawgraph(nls,yrange,method_list,'med_r','Median Rotation Error',...
    'Number of Points','Rotation Error (degrees)');

figure('color','w','position',[w*i,100,w,h]);i=i+1;
subplot(121)
xdrawgraph(nls,yrange,method_list,'mean_t','Mean Translation Error',...
    'Number of Points','Translation Error (%)');
subplot(122)
xdrawgraph(nls,yrange,method_list,'med_t','Median Translation Error',...
    'Number of Points','Translation Error (%)');


function xdrawgraph(xs,yrange,method_list,field,ti,lx,ly)

box('on');
hold('all');

p= zeros(size(method_list));
for i= 1:length(method_list)
    p(i)= plot(xs,method_list(i).(field),method_list(i).marker,...
        'color',method_list(i).color,...
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

end
