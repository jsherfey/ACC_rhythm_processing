% Model: TC-RE

spec=[];
spec.nodes(1).label = 'TC';
spec.nodes(1).multiplicity = 1;
spec.nodes(1).dynamics = {'V''=current/c'};
spec.nodes(1).mechanisms = {'Ca_TC','iH_TC','iK_TC','iKLeak_TC','iLeak_TC','iNa_TC','iT_TC','itonic'};
spec.nodes(1).parameters = {'c',1,'V_IC',-65};
spec.nodes(2).label = 'RE';
spec.nodes(2).multiplicity = 1;
spec.nodes(2).dynamics = {'V''=current/c'};
spec.nodes(2).mechanisms = {'iK_RE','iLeak_RE','iNa_RE','iT_RE','itonic'};
spec.nodes(2).parameters = {'c',1,'v_IC',-85};
spec.connections(1,2).label = 'TC-RE';
spec.connections(1,2).mechanisms = {'iAMPA'};
spec.connections(1,2).parameters = {'gAMPA',0.69};
spec.connections(2,1).label = 'RE-TC';
spec.connections(2,1).mechanisms = {'iGABAA','iGABAB'};
spec.connections(2,1).parameters = {'gGABAA_base',0.069,'gGABAB',0.138}
spec.connections(2,2).label = 'RE-RE';
spec.connections(2,2).mechanisms = {'iGABAA'};
% spec.connections(2,2).mechanisms = {'iAMPA_PY_faux','iGABAA'};
spec.connections(2,2).parameters = {'gGABAA_base',0.138};
% spec.connections(2,2).parameters = {'gAMPA',0.1,'PY_square_freq',0.0,'gGABAA_base',0.138};
% spec.connections(1,1).label = 'TC-TC';
% spec.connections(1,1).mechanisms = {'iAMPA_PY_faux'};
% spec.connections(1,1).parameters = {'gAMPA',0.1,'PY_square_freq',0.0};

codepath='~/x010-dnsim/dnsim';
plotvars_flag = 1;
plot_flag = 1;
plotpower_flag = 1;
plotpacoupling_flag = 1;
saveplot_flag = 1;
overwrite_flag = 1;
timelimits = [0 5000];
SOLVER = 'euler';
rootdir = 'DATA DIRECTORY GOES HERE'
addpath = codepath;
cluster_flag = 0;
% cluster_flag = 1;
savedata_flag = 1;

[specs,timestamp,rootoutdir]=simstudy(spec,{'(TC,RE)','TC'},{'stim','gH'},{'[0,-0.1,-0.3,-0.5,-0.7,-1.0,-1.5,-2.0,-2.5]','[1.0,0.5,0.3,0.1,0.07,0.05,0.01,0.007,0.005,0.0005]'},'plotvars_flag',plotvars_flag,'plot_flag',plot_flag,'plotpower_flag',plotpower_flag,'plotpacoupling_flag',plotpacoupling_flag,'saveplot_flag',saveplot_flag,'overwrite_flag',overwrite_flag,'timelimits',timelimits,'SOLVER',SOLVER,'rootdir',rootdir,'addpath',addpath,'cluster_flag',cluster_flag,'savedata_flag',savedata_flag);
% [specs,timestamp,rootoutdir]=simstudy(spec,{'(TC,RE)','(TC,RE)','(TC,RE)','TC','RE-TC'},{'stim','cos_ampl','cos_freq','gH','spm'},{'[0, -0.1, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0]','[0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]','[0.6, 1.0, 1.3]','[1.0, 0.5, 0.1, 0.07, 0.05, 0.01, 0.007, 0.005, 0.0005]','[2,3]'},'plotvars_flag',plotvars_flag,'plot_flag',plot_flag,'plotpower_flag',plotpower_flag,'plotpacoupling_flag',plotpacoupling_flag,'saveplot_flag',saveplot_flag,'overwrite_flag',overwrite_flag,'timelimits',timelimits,'SOLVER',SOLVER,'rootdir',rootdir,'addpath',addpath,'cluster_flag',cluster_flag,'savedata_flag',savedata_flag);

exit
