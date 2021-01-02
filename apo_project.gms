*title APO Refrigerant Optimisation

Sets
 i row labels /CH3, CH2, CH, C, CH2DCH, CHDCH, CH2DC, CHDC, CDC, CH2DCDCH, CH2DCDC, CHDCDCH, CH3O, CH2O, CHDO, CDO, CH2NH2, CHNH2, CNH2, CH3NH, CH2NH, CHNH, CH3N, CH2N, CH2F, CHF, CF, CHF2, CF2, CF3/
 j column labels /Tmi, Tbi, Tci, Pci, Hvi, Cpi, Val, Nmax/
 k 'number index' /1*3/;
*removed set m and included equations for 272 and 316 for all equations using m
Parameters
 Tm0 'Melting point constant [K]' /147.45/
 Tb0 'Boiling point constant [K]' /222.543/
 Tc0 'critical temperature constant [K]' /231.239/
 Pc01 'critical pressure constant 1 [bar]' /5.9827/
 Pc02 'critical pressure constant 2 [bar^-0.5]' /0.108998/
 Hv0 'enthalpy of vaporisation constant [kJ/mol]' /11.733/
 HvR134a 'heat of vaporisation for R134a at 272 K [kJ/mol]' /20.33/
 CpR134a 'heat capacity of R134a at 294 K [kJ/kmol K]' /143.9/

*Were these meant to be 10, 21, 32, 43? Ask Sochi
*Don't think you can define exponents like this - throws error.
 x(k) 'exponents'
 / 1 0
   2 1
   3 2 /

 thermodynamics(i,j);

$onUNDF
$Call GDXXRW apo_project.xlsx trace=3 par=thermodynamics rng=Sheet2!A1:H31 rdim=1 cdim=1
$GDXIN apo_project.gdx
$LOAD  thermodynamics
$GDXIN

Display thermodynamics;

Positive Variables
 Tm 'Melting point [K]'
 Tb 'Boiling pouint [K]'
 Tc 'Critical temperature [K]'
 Hv298 'heat of vaporisation at 298 K [kJ/mol]'
 Hv272 'heat of vaporisation at 272 K [kJ/mol]'
 Pc 'Critical Pressure [bar]'
 Cp 'Heat Capacity [kJ/kmol K]'
 Tbr 'Temperature ratio Tb/Tc'
 Tr_272 'Temperature ratio at 276 K'
 Tr_316 'Temperature ratio at 316 K'
 Pv_272 'Vapour pressure at 272 K [bar]'
 Pv_316 'Vapour pressure at 316 K [bar]'

Variables
 W 'acentric factor'
 f0Tbr 'Pitzer term 1 with Tr=Tbr'
 f1Tbr 'Pitzer term 2 with Tr=Tbr'

 f0_272  'Pitzer term 1 at 272 K'
 f1_272  'Pitzer term 2 at 272 K'
 f2_272  'Pitzer term 3 at 272 K'

 f0_316  'Pitzer term 1 at 316 K'
 f1_316  'Pitzer term 2 at 316 K'
 f2_316  'Pitzer term 3 at 316 K'
 z  'objective variable';

*Default lower bound for integer variables is already 0.
Integer variable
 N(i) 'Number of group i';


Binary variables
 y(k,i) 'binary';

Equations
 Num(i) 'Number of each group'
 Tmelt  'Melting point'
 Tboil  'Boiling point'
 Tcrit  'Critical temperature'
 Hvap298 'heat of vaporisation at 298 K'
 Cpsum  'heat capacity'
 Pcsum  'Critical pressure'

 Hvap272 'heat of vaporisation at 272 K'

 Tbrat  'Boiling temperature ratio'
 Trat_272 'Temperature ratio at 276 K'
 Trat_316 'Temperature ratio at 316 K'

 f0AW_272 'Ambrose Walton equation for f0 at 272 K'
 f1AW_272 'Ambrose Walton equation for f1 at 272 K'
 f2AW_272 'Ambrose Walton equation for f2 at 272 K'

 f0AW_316 'Ambrose Walton equation for f0 at 316 K'
 f1AW_316 'Ambrose Walton equation for f1 at 316 K'
 f2AW_316 'Ambrose Walton equation for f2 at 316 K'

 f0TbAW  'Ambrose Walton equation for f0 at the boiling point'
 f1TbAW  'Ambrose Walton equation for f1 at the boiling point'
 Ace  'Acentric factor'

 Pvcorr_272  'Vapour pressure correlation at 272 K'
 Pvcorr_316  'Vapour pressure correlation at 316 K'

*objective
 obj 'objective function'

*Process/physical constraints
 Pv_con272 'Vapour pressure constraint at 272 K'
 Pv_con316 'Vapour pressure constraint at 316 K'

*Pvc 'Vapour pressure constraint at 316 K'
 enthalpy 'enthalpy of vaporisation constraint'
 heatcap  'heat capacity constraint'
 phase  'no-solid constraint'

*Structural constraints
 totalgroups  'constraint on total number of groups'
 maxgroup(i) 'limits number of each group'
 valency  'check on valency'
 minbonds 'minimum bonds allowed'
 maxbonds 'max bonds allowed';
*nextjoin 'disallows adjacent group double bonding';

*integer cuts
*number of each group as binary combination
 Num(i).. N(i)=e= sum(k, (2**x(k))*y(k,i));

*Calculating Parameters
 Tmelt..  Tm =e= Tm0*log(sum(i, N(i)*thermodynamics(i, 'Tmi')));
 Tboil..  Tb =e= Tb0*log(sum(i, N(i)*thermodynamics(i, 'Tbi')));
 Tcrit..  Tc =e= Tc0*log(sum(i, N(i)*thermodynamics(i, 'Tci')));
 Hvap298.. Hv298 =e= Hv0 + sum(i, N(i)*thermodynamics(i, 'Hvi'));
 Cpsum.. Cp =e= sum(i, N(i)*thermodynamics(i, 'Cpi'));
 Pcsum.. Pc =e= Pc01 +(Pc02 + sum(i, N(i)*thermodynamics(i, 'Pci')))**(-2);

*Adjusting Parameters
 Hvap272.. Hv272 =e= Hv298*(((Tc-272)/(Tc-298))**0.38);

*Defining temperature ratios - Tr(272) and Tr(316) not recognised
 Tbrat..       Tbr =e= Tb/Tc;
 Trat_272.. Tr_272 =e= 272/Tc;
 Trat_316.. Tr_316 =e= 316/Tc;

*Pitzer expansion parameters at 272K
 f0AW_272.. f0_272 =e= (-5.97616*(1-Tr_272)+1.29874*(1-Tr_272)**1.5 -0.60394*(1-Tr_272)**2.5 -1.06841*(1-Tr_272)**5)/ Tr_272;
 f1AW_272.. f1_272 =e= (-5.03365*(1-Tr_272)+1.11505*(1-Tr_272)**1.5 -5.41217*(1-Tr_272)**2.5 -7.46628*(1-Tr_272)**5)/ Tr_272;
 f2AW_272.. f2_272 =e= (-0.64771*(1-Tr_272)+2.41539*(1-Tr_272)**1.5 -4.26979*(1-Tr_272)**2.5 +3.25259*(1-Tr_272)**5)/ Tr_272;

*Pitzer expansion parameters at 316K
 f0AW_316.. f0_316 =e= (-5.97616*(1-Tr_316)+1.29874*(1-Tr_316)**1.5 -0.60394*(1-Tr_316)**2.5 -1.06841*(1-Tr_316)**5)/ Tr_316;
 f1AW_316.. f1_316 =e= (-5.03365*(1-Tr_316)+1.11505*(1-Tr_316)**1.5 -5.41217*(1-Tr_316)**2.5 -7.46628*(1-Tr_316)**5)/ Tr_316;
 f2AW_316.. f2_316 =e= (-0.64771*(1-Tr_316)+2.41539*(1-Tr_316)**1.5 -4.26979*(1-Tr_316)**2.5 +3.25259*(1-Tr_316)**5)/ Tr_316;

*Acentric factor
 f0TbAW.. f0Tbr =e= (-5.97616*(1-Tbr)+1.29874*(1-Tbr)**1.5 -0.60394*(1-Tbr)**2.5 -1.06841*(1-Tbr)**5)/ Tbr;
 f1TbAW.. f1Tbr =e= (-5.03365*(1-Tbr)+1.11505*(1-Tbr)**1.5 -5.41217*(1-Tbr)**2.5 -7.46628*(1-Tbr)**5)/ Tbr;
 Ace.. W =e= -(log(Pc/1.01325)+ f0Tbr)/f1Tbr ;

*Vapour Pressure
 Pvcorr_272.. Pv_272 =e= Pc* exp(f0_272+ W*f1_272 + (W**2)*f2_272);
 Pvcorr_316.. Pv_316 =e= Pc* exp(f0_316+ W*f1_316 + (W**2)*f2_316);

 obj.. z =e= Cp/Hv272;

*Physical constraints - PV(272) and PV(316) are NOT defined
 Pv_con272.. Pv_272 =G= 1.1;
 Pv_con316.. Pv_316 =L= 14;
 enthalpy.. Hv272 =G= 20.33;
 heatcap.. Cp =L= 143.9;
 phase.. Tm =L= 272;

*Constraints on molecule formation
 totalgroups.. sum(i, N(i))=L= 15;
 maxgroup(i).. N(i) =L= 5;


 valency..  sum(i,(2-thermodynamics(i,'Val'))*N(i))=e= 2;
 minbonds.. sum(i,N(i)*thermodynamics(i,'Val'))=G= 2*(sum(i,N(i))-1);
 maxbonds.. sum(i,N(i)*thermodynamics(i,'Val'))=L= (sum(i, N(i)))*(sum(i,N(i))-1);
*nextjoin.. N(i)*(thermodynamics(i, 'Val')-1) + 2 - sum(i, N(i)) =L= 0;

*Solving
Model molecule /all/;
Solve molecule using minlp minimizing z;
