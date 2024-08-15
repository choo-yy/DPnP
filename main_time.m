clear; clc;
addpath epnp;
addpath lhm;
addpath func;
addpath dls_pnp_matlab;
addpath RDLT;
addpath OPnP;
addpath DPnP;

% experimental parameters
ntest= 100;
nl= 5;
npts= 10:40:1040;

name=   {'DLT', 'EPnP+GN', 'RDLT+GN', 'DLS', 'LHM', 'RPnP',     'OPnP', 'DPnP', 'DPnP+GN'};
f=      { @DLT,  @EPnP_GN,  @RDLT_GN,  @DLS,  @LHM,  @RPnP,  @OPnP1res,  @DPnP,  @DPnP_GN};
marker= { '-*',      '-o',      '-x',  '->',  '-+',   '-^',       '-s',   '-d',      '-p'};


method_list= struct('name', name, 'f', f, 't', zeros(size(npts)),...
    'marker', marker);

% camera's parameters
width= 640;
height= 480;
f= 800;

for j= 1:length(npts)
    npt= npts(j);
    fprintf('%d points\n',npt);
    
    % generate experimental data
    para= cell(ntest,2);
    for i= 1:ntest
        % generate 3d coordinates in target space
        XYZrng= [-0.5 0.5]*height;
        XXw= xrand(3,npt,XYZrng);
        % generate rotation matrix and translation vector
        R= rodrigues(randn(3,1));
        t= [(rand-0.5)*width; (rand-0.5)*height; f]*rand*2 + [0; 0; f];
        % project to image points
        XX= R*XXw+repmat(t,1,size(XXw,2));
        xx= [XX(1,:)./XX(3,:); XX(2,:)./XX(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;
		xxn= xxn/f;
        % save
        para{i,1}= XXw;
        para{i,2}= xxn;
    end

    for k= 1:length(method_list)
        tic;
        for i= 1:ntest
            XXw= para{i,1};
            xxn= para{i,2};
            method_list(k).f(XXw,xxn);
        end
        t= toc; method_list(k).t(j)= t/ntest*1000;
        disp([method_list(k).name ' - ' num2str(t) ' s']);
    end


end
%%
close all;

set(groot,'defaultLineLineWidth',2)
% set(groot,"DefaultAxesFontSize",18)
figure('color','w','position',[100 100 360 320]);
hold all;
box on;
p= zeros(size(method_list));
for k= 1:length(method_list)
    p(k)= plot(npts,method_list(k).t,...
        method_list(k).marker,...
        'displayname',method_list(k).name);
end
legend(p);
xlim(npts([1,end]));
ylim([0 10]);
% xtick= 4:12:100;
% set(gca,'xtick',xtick);

xlabel('Number of Points');
ylabel('Computational Time (milliseconds)');

