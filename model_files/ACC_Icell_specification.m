function pop = ACC_Icell_specification
% Generic ACC I-cell model (CB+ or PV+ interneuron)

pop.equations='dv/dt=@current+input; input=0; {iNa,iK}; v(0)=-65';
pop.parameters={};
pop.name='';
pop.size=2;
pop.mechanism_list={};
