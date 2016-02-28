function model = ACC_Icell_model(varargin)
% Generic ACC I-cell model (CB+ or PV+ interneuron)

s.pops.equations='dv/dt=@current+input; input=0; {iNa,iK}';
if nargin>0
  s.pops.parameters=varargin;
end

model=GenerateModel(s,'open_link_flag',1);
