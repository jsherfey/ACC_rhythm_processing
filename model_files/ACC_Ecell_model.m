function model = ACC_Ecell_model(varargin)
% ACC E-cell model (Pyramidal cell)

s.pops.equations='dv/dt=@current+input; input=0; {iNa,iK}';
if nargin>0
  s.pops.parameters=varargin;
end

model=GenerateModel(s,'open_link_flag',1);
