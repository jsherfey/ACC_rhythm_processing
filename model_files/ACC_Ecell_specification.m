function pop = ACC_Ecell_specification
% ACC E-cell model (Pyramidal cell)

pop.equations='dv/dt=@current+input; input=0; {iNa,iK}; v(0)=-65';
pop.parameters={};
pop.name='';
pop.size=2;
pop.mechanism_list={};
