clear; clc;
addpath epnp;
addpath lhm;
addpath func;
addpath dls_pnp_matlab;
addpath RDLT;
addpath OPnP;
addpath DPnP;

% experimental parameters
nl= 2;
npts= 4:20;
num= 1000;

% compared methods
A= zeros(size(npts));
B= zeros(num,1);
name=   {'DLT', 'EPnP+GN', 'RDLT+GN', 'DLS', 'LHM', 'RPnP',     'OPnP', 'DPnP', 'DPnP+GN'};
f=      { @DLT,  @EPnP_GN,  @RDLT_GN,  @DLS,  @LHM,  @RPnP,  @OPnP1res,  @DPnP,  @DPnP_GN};
marker= { '-*',      '-o',      '-x',  '->',  '-+',   '-^',       '-s',   '-d',      '-p'};


method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker);

% experiments
for i= 1:length(npts)
    
    npt= npts(i);
    fprintf('npt = %d: ',npt);
    
    for j= 1:num
        
        % camera's parameters
        width= 640;
        height= 480;
        f= 800;
        
        % generate 3d coordinates in camera space
        % Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])]; % ordinary 3d
        Xc= [xrand(1,npt,[1 2]); xrand(1,npt,[1 2]); xrand(1,npt,[4 8])]; % quasi-singular
        t= mean(Xc,2) + rand(3,1);
        R= rodrigues(randn(3,1));
        XXw= inv(R)*(Xc-repmat(t,1,npt));
        
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
yrange= [0 3];

i= 0; w= 1000; h= 350;
set(groot,'defaultLineLineWidth',2)
figure('color','w','position',[w*i,100,w,h]);%i=i+1;
subplot(121)
xdrawgraph(npts,yrange,method_list,'mean_r','Mean Rotation Error',...
    'Number of Points','Rotation Error (degrees)');
subplot(122)
xdrawgraph(npts,yrange,method_list,'med_r','Median Rotation Error',...
    'Number of Points','Rotation Error (degrees)');

figure('color','w','position',[w*i,100,w,h]);%i=i+1;
subplot(121)
xdrawgraph(npts,yrange,method_list,'mean_t','Mean Translation Error',...
    'Number of Points','Translation Error (%)');
subplot(122)
xdrawgraph(npts,yrange,method_list,'med_t','Median Translation Error',...
    'Number of Points','Translation Error (%)');

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

end