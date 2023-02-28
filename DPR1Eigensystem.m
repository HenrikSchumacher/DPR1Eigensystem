(* ::Package:: *)

BeginPackage["DPR1Eigensystem`",{"CCompilerDriver`"}]

ClearAll[DPR1Apply];
DPR1Apply::usage="DPR1Apply[diag,z,x] computes the action of the matrix DiagonalMatrix[diag]+KroneckerProduct[z,z] on the vector of list of vectors x  in a more efficient way.";

ClearAll[DPR1TestEigensystem];
DPR1TestEigensystem::usage="DPR1Apply[diag,z,\[Lambda],u] computes relative error of Max[Abs[(DiagonalMatrix[diag]+KroneckerProduct[z,z])/u - \[Lambda] u)]]/\[Lambda].";

ClearAll[DPR1ApplyShiftedInverse];
DPR1ApplyShiftedInverse::usage="DPR1ApplyShiftedInverse[diag,z,\[Mu],x,nf] computes LinearSolve[DiagonalMatrix[diag-\[Mu]]+KroneckerProduct[z,z],x] in a more efficient way. Use nf = Nearest[diag->\"Index\"].";


ClearAll[DPR1RaleighIterations];
DPR1RaleighIterations::usage="DPR1RaleighIterations[diag,z,\[Mu],x,nf] attempts to employ Raleigh iteration with shift \[Mu] to find an eigenvalue of the matrix DiagonalMatrix[diag]+KroneckerProduct[z,z] that lies close to \[Mu]. Use nf = Nearest[diag->\"Index\"].";

ClearAll[DPR1EigenvectorFromEigenvalue];
DPR1EigenvectorFromEigenvalue::usage="DPR1EigenvectorFromEigenvalue[diag,z,\[Mu]] computes the eigenvector(s) of the matrix DiagonalMatrix[diag]+KroneckerProduct[z,z] corresponding to the eigenvalue(s) \[Mu].";

ClearAll[DPR1Eigenvalues];
DPR1Eigenvalues::usage="DPR1Eigenvalues[diag,z] computes the eigenvalues of the matrix DiagonalMatrix[diag]+KroneckerProduct[z,z].";

ClearAll[DPR1Eigensystem];
DPR1Eigensystem::usage="DPR1Eigensystem[diag,z] computes the eigenvalues as well as the eigenvectors of the matrix DiagonalMatrix[diag]+KroneckerProduct[z,z].";

Begin["`Private`"];

DPR1Apply[diag_?VectorQ,z_?VectorQ,x_]:=cDPR1Apply[diag,z,x];

cDPR1Apply=Compile[{{diag,_Real,1},{z,_Real,1},{x,_Real,1}},
	Block[{n,sum,xi,y},
		n=Length[diag];
		sum=0.;
		y=Table[
			xi=Compile`GetElement[x,i];
			sum+=Compile`GetElement[z,i]xi;
			Compile`GetElement[diag,i] xi
		,{i,1,n}];
		
		Do[y[[i]]+=Compile`GetElement[z,i]sum;,{i,1,n}];
		y
	],
	CompilationTarget->"C",
	RuntimeAttributes->{Listable},
	Parallelization->True,
	RuntimeOptions->"Speed"
];

DPR1TestEigensystem[diag_?VectorQ,z_?VectorQ,\[Lambda]_,U_]:=cDPR1TestEigensystem[diag,z,\[Lambda],U];

cDPR1TestEigensystem=Compile[{{diag,_Real,1},{z,_Real,1},{\[Lambda],_Real},{u,_Real,1}},
	Block[{n,sum,ui,v,max,maxerror},	
		maxerror=0.;
		n=Length[diag];
		sum=0.;
		v=Table[
			ui=Compile`GetElement[u,i];
			sum+=Compile`GetElement[z,i]ui;
			Compile`GetElement[diag,i] ui
		,{i,1,n}];
		
		Do[
			maxerror=Max[
			maxerror,
			Abs[Compile`GetElement[v,i]+Compile`GetElement[z,i]sum-\[Lambda] Compile`GetElement[u,i]]
		];
		,{i,1,n}];
		
		maxerror/Abs[\[Lambda]]
	],
	CompilationTarget->"C",
	RuntimeAttributes->{Listable},
	Parallelization->True,
	RuntimeOptions->"Speed"
];

DPR1ApplyShiftedInverseNoResonance[diag_?VectorQ,z_?VectorQ,\[Mu]_,x_?VectorQ]:=Module[{y,yz,OnePluszyz},
	y=1./(diag-\[Mu]);
	yz=y z;
	OnePluszyz=1.+z . yz;
	If[Abs[OnePluszyz]>10. $MachineEpsilon,
		y x-yz(yz . x/OnePluszyz)
		,
		(*Not sure what to do here. At least just returning x seems to work for DPR1RaleighIterations.*)
		x
	]
];

ClearAll[DPR1ApplyShiftedInverseResonance];
DPR1ApplyShiftedInverseResonance[diag_?VectorQ,z_?VectorQ,k_,x_?VectorQ]:=Module[{y,yz,w,b,tmp,u},
	y=Join[1./(diag[[;;k-1]]-diag[[k]]),{0.},1./(diag[[k+1;;]]-diag[[k]])];
	yz=z y;
	w=yz/(-z[[k]]);
	b=(1.+z . yz)/z[[k]]^2;
	tmp=w . x+b x[[k]];
	u=y x+w x[[k]];
	u[[k]]+=tmp;
	u
];

DPR1ApplyShiftedInverse[diag_?VectorQ,z_?VectorQ,\[Mu]_,x_?VectorQ,nf_NearestFunction,TOL_:10. $MachineEpsilon]:=Module[{k},
	(*Applies the inverse of the matrix DiagonalMatrix[diag-\[Mu]]+KroneckerProduct[z,z] to the vector x.*)
	
	(*We need to detect whether the shifted diagonal diag-\[Mu] 
	has a zero entry. We use a NearestFunction to speed this up.*)
	
	(*Find the entry closest to the shift \[Mu].*)
	k=nf[\[Mu]][[1]];
	
	If[
		Abs[diag[[k]]-\[Mu]]>TOL Abs[\[Mu]]
		,
		DPR1ApplyShiftedInverseNoResonance[diag,z,\[Mu],x]
		,
		DPR1ApplyShiftedInverseResonance[diag,z,k,x]
	]
];

ClearAll[DPR1RaleighIterations];
DPR1RaleighIterations[diag_?VectorQ,z_?VectorQ,\[Mu]0_,nf_NearestFunction]:=Module[{n,\[Mu],\[Mu]old,iter,u,v,s,\[Epsilon],maxiter},
	\[Epsilon] = 100 $MachineEpsilon;
	maxiter = 1000;
	
	n=Length[diag];
	
	\[Mu]=\[Mu]0;
	\[Mu]old=2.\[Mu];
	iter=0;
	
	(*Compute a rough initial guess for the eigenvector*)
	u=Table[
		s=diag[[k]]-\[Mu];
		If[Abs[s]>\[Epsilon] Abs[\[Mu]],z[[k]]/s,0.]
	,{k,1,n}];
	
	
	While[(Abs[\[Mu]-\[Mu]old]>Abs[\[Mu]]\[Epsilon])&&(iter<maxiter),
		++iter;
		
		\[Mu]old=\[Mu];
		
		(* u = LinearSolve[A - \[Mu] IdentityMatrix[n],u ] *)
		
		u=Normalize[DPR1ApplyShiftedInverse[diag,z,\[Mu],u,nf]];
		
		(*Raleigh quotient \[Mu] = u.A.u / Norm[u] *)
		\[Mu]=(u . DPR1Apply[diag,z,u]);
	];
	\[Mu]
	
];

DPR1EigenvectorFromEigenvalue[diag_?VectorQ,z_?VectorQ,\[Lambda]_] := cDPR1EigenvectorFromEigenvalue[diag,z,\[Lambda]];

cDPR1EigenvectorFromEigenvalue = Compile[{{diag,_Real,1},{z,_Real,1},{\[Lambda],_Real}},
	Block[{n,u,x,sum,invnorm},
		n=Length[diag];
		sum=0.;
		u=Table[
			x=Compile`GetElement[z,i]/(Compile`GetElement[diag,i]-\[Lambda]);
			sum+= x x;
			x
		,{i,1,n}];
		
		invnorm=1./Sqrt[sum];
		Do[u[[i]]=invnorm Compile`GetElement[u,i];,{i,1,n}];
		u
	],
	CompilationTarget->"C",
	RuntimeAttributes->{Listable},
	Parallelization->True,
	RuntimeOptions->"Speed"
];

Options[DPR1Eigenvalues]={
	"Tolerance"->100 $MachineEpsilon,
	"MaxIterations"->20
};

DPR1Eigenvalues[diag_?VectorQ,z_?VectorQ,OptionsPattern[]]:=Module[{p,zz,a0,b0},
	p=Ordering[diag];
	zz=(z z)[[p]];
	a0=diag[[p]];
	b0=Append[Rest[a0],a0[[-1]]+ Total[zz]];
	Reverse[cDPR1Eigenvalues[a0,zz,a0,b0,OptionValue["MaxIterations"],OptionValue["Tolerance"]]]
];

cDPR1Eigenvalues=With[{
	eps=N[$MachineEpsilon]
},
Compile[{{diag,_Real,1},{zz,_Real,1},{a0,_Real},{b0,_Real},{maxiter,_Integer},{RelTOL,_Real}},
	Module[{n,a,fa,b,fb,c,fc,Dfc,TOL,cold,fcold,x,\[Delta]c,iter,\[Epsilon],sum,\[Delta]},
		n=Length[diag];
		a=a0;
		b=b0;
		\[Epsilon]=10.eps;
		(*TOL=RelTOL(b-a);*)
		
		TOL=RelTOL Max[Abs[b],Abs[a]];
		
		c=0.5(a+b);
		
		(*Evaluate the characteristic polynomial f at c.*)
		sum=0.;
		Do[
			\[Delta]=Compile`GetElement[diag,i]-c;
			sum+=Compile`GetElement[zz,i]/\[Delta];
		,{i,1,n}];
		fc=1.+sum;
		
		cold=c;
		fcold=fc;
		
		(*First we have to bracket the solution.*)
		(*Trying to find {a,b} with f[a] < 0 < f[b].*)
		If[ fc >= 0.,
		(
			While[fc>=0.,
			cold=c;fcold=fc;
			c=0.5(a +c);
			
			(*Evaluate the characteristic polynomial f at c.*)
			sum=0.;
			Do[
				\[Delta]=Compile`GetElement[diag,i]-c;
				sum+=Compile`GetElement[zz,i]/\[Delta];
			,{i,1,n}];
			fc=1.+sum;
			];
			
			a=c;fa=fc;
			b=cold;fb=fcold;
		),(
			While[fc<0.,
			cold=c;fcold=fc;
			c=0.5(c+b);
		
			(*Evaluate the characteristic polynomial f at c.*)
			sum=0.;
			Do[
				\[Delta]=Compile`GetElement[diag,i]-c;
				sum+=Compile`GetElement[zz,i]/\[Delta];
			,{i,1,n}];
			fc=1.+sum;
			];
			
			a=cold;fa=fcold;
			b=c;fb=fc;
		)];
		
		(*Refinement.*)
		c=0.5(a+b);
		
		(*Evaluate the characteristic polynomial f at c.*)
		sum=0.;
		Do[
			\[Delta]=Compile`GetElement[diag,i]-c;
			sum+=Compile`GetElement[zz,i]/\[Delta];
		,{i,1,n}];
		fc=1.+sum;
		
		iter=0;
		
		While[iter<maxiter,
			++iter;
			
			(*First try a Newton step.*)
			(*\[Delta]c=-fc/cDf[diag,zz,c];*)
			
			(*Evaluate the derivative Dfc of the characteristic polynomial f at c.*)
			Dfc=0.;
			Do[
				\[Delta]=Compile`GetElement[diag,i]-c;
				Dfc+=Compile`GetElement[zz,i]/(\[Delta] \[Delta]);
			,{i,1,n}];
			\[Delta]c=-fc/Dfc;
			
			x=c+\[Delta]c;
			If[a<x<b,
			(
				(*Newton's method gives a valid iterate.*)
				If[Abs[\[Delta]c]<Max[\[Epsilon],TOL],
				(*Print["A"->a<x<b];*)
				Return[c];
				];
				c=x;
				
				(*Evaluate the characteristic polynomial f at c.*)
				sum=0.;
				Do[
					\[Delta]=Compile`GetElement[diag,i]-c;
					sum+=Compile`GetElement[zz,i]/\[Delta];
				,{i,1,n}];
				fc=1.+sum;	
			),(
				(*Newton's method leaves the bracketing interval; we continue with the interval's midpoint.*)
				If[fc<0.,
				a=c;fa=fc;
				,
				b=c;fb=fc;
				];
				
				If[Abs[b-a]<Max[\[Epsilon],TOL], Return[c];];
				
				c=0.5(a+b);
				
				(*Evaluate the characteristic polynomial f at c.*)
				sum=0.;
				Do[
					\[Delta]=Compile`GetElement[diag,i]-c;
					sum+=Compile`GetElement[zz,i]/\[Delta];
				,{i,1,n}];
				fc=1.+sum;
			)];
		];
		
		(*We ran out of iterations. Return the result anyways.*)
		
		Return[c]
	],
	CompilationTarget->"C",
	RuntimeAttributes->{Listable},
	Parallelization->True,
	RuntimeOptions->"Speed"
]
];

Options[DPR1Eigensystem]={
	"Tolerance" -> 100 $MachineEpsilon,
	"MaxIterations" -> 20
};

DPR1Eigensystem[diag_?VectorQ,z_?VectorQ,OptionsPattern[]]:=Module[{\[Lambda],U},
	\[Lambda]=DPR1Eigenvalues[diag,z,
		"MaxIterations"->OptionValue["MaxIterations"],"Tolerance"->OptionValue["Tolerance"]
	];
	
	U=DPR1EigenvectorFromEigenvalue[diag,z,\[Lambda]];
	
	{\[Lambda],U}
];

End[];

EndPackage[];
