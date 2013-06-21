(* Mathematica Package *)
(* Created by the Wolfram Workbench Sep 3, 2010 *)
(* This package contains core and support functions to determine the entropy of the outputs of a Binary Partition Network (see Billings et al)
Note that efficient computation should make use to the assiciated C library. Mathematca functions for computing entropy are also included however*)

BeginPackage["BinaryPartitionNetworks`"]

(*CORE FUNCTIONS********************************************************************)

p\[Phi]::usage = "p\[Phi][\[Mu],\[Beta],\[Phi],d]: Gives the probability that an output is active if that output makes d random connections to \[Mu] inputs \ 
where \[Beta] of those inputs are active and the output has a threshold of \[Phi]."

Pout::usage = "Pout[d, pin, \[Phi]]: Proababilty that an output is active in a BPN where outputs make d connections to inputs and have a thereshold \
\[Phi]. Inputs active with probability pin."

Joint\[Alpha]\[Beta]::usage = "Joint\[Alpha]\[Beta][\[Mu], pon, \[Gamma], d, \[Phi]]: Gives the joint probabilty of output success class \[Alpha] and \
input success class \[Beta], where there are \[Mu] inputs active with probability pon and each output makes d connections with threshold \[Phi]."

Uniquecodes::usage = "Uniquecodes[indist,\[Mu],tc]: Gives a vector of length \[Mu]+1 containing the number of unique input codes in each of the \[Mu]+1 \
possible input success classes with associated p.d.f indist given a sample of tc events (input realisations)."

Binmat::usage="Binmat[done, dtwo, s]: Random binary adjacency matrix for a BPN having done outputs, dtwo inputs and s connections per output."

connectedMscheme::usage = "connectedMscheme[\[Mu],d,\[Delta],\[Epsilon],mdims]: Generates a transition matrix for an ensemble of cells who connect to their \
inputs by drawing 1 input from \[Mu] given that \[Delta] have already been drawn and replacement is not permitted. mdims must be \[Mu]-d+1, and is the length\
of the support for the connected component p.d.f. \[Epsilon] permits analysis of the matrix but should normally be 0. "

Occupancies::usage = "Occupancies[nmat,ini,\[Gamma]]: Returns the connected component p.d.f for each output up to a total of \[Gamma] outputs. nmat is the \
appropriate transition matrix found with 'connectedMscheme' and ini is a vector of length \[Mu]-d+1, containing the initial connected component \
p.d.f."

EntropyApprx::usage = "Entropyapprx[tc,\[Mu],d,\[Phi],pon,\[Gamma]]: Approximate entropy of the outputs of a Binary Partition Network in which \[Gamma] \
outputs make d connections per output with \[Mu] inputs active with probability pon, where each output has threshold \[Phi] and tc events are encoded \
by the network."

EntropyapprxDEP::usage="Deprecated version of Entropyapprx. Is very slow."

Bitent::usage = "Bitent[p]: Information entropy of a single Bernoulli trial, succesful with probability p."

Sourceent::usage = "Sourceent[\[Mu],pon,uc]: Source entropy of \[Mu] independent units active with probability pon, given that uc[[x]] unique codes occur in input class x.\
when events are encoded. uc must be a vector of length \[Mu]+1."

Typset::usage = "Typset[r,H,\[Mu]] gives the number of patterns within a fraction r of the typical set associated with \[Mu] random inputs active with probability p."

Levycost::usage = "Levycost[\[Mu],pin,x,\[Gamma],pout,y,\[Delta]]: Basic model of the running costs of a BPN having \[Mu] inputs active with probability pin and \[Gamma] \
outputs active with probability pout, where the ratio of energy expenditure of active to inactive inputs is x and the ratio of expenditure of \
active to inactive outputs is y.  \[Delta] is the ratio of the cost of an off output to an off input. Measure is in units of the energy expended by an off input. \
Based on model due to Levy & Baxter, 1996, Neural Computation Vol. 8, No. 3, Pages 531-543. Cost = \[Mu][pin*x+(1-pin)]+\[Delta]*\[Gamma][pout*y+(1-pout)]."

Levycost2::usage = "Levycost2[\[Mu],pin,x,\[Gamma],pout,y,\[Delta],G,d,connscal]: This model is as Levycost, but includes a term that accounts for the cost associated\
with the activity of the glomeruli"

Hdist::usage="Hdist[\[Mu], d, \[Phi], \[Gamma], pon]: Exact determination of (hout-hin)/\[Mu]*pon where hout and hin are the expected Hamming distances between \ 
pairs of output and input binary vectors respectively. mode= 1. Random binary vectors, 2. Network outputs, 3. Random vectors with same sparsity as \
network outputs" 

HdistNUM::usage="HdistNUM[\[Mu], d, \[Phi], \[Gamma], pon, pairs]: Numerical determination of (hout-hin)/\[Mu]*pon where hout and hin are the expected Hamming distances between\
pairs of output and input binary vectors respectively." 

GThreshold::usage="Threshold[pin_, d_, gain_, offset_: 0] Maps intervals of pin into threshold values. Offset=0: No tonic inhibition. >0, Floor[d/2] offset."

Thrent::usage="Thrent[data_, pmin_, pinc_, pmax_, d_, gain_, offset_: 0] Entropy encoded when threshold increases with increasing pin."

Thrprune::usage="Thrprune[data_,pmin_, pinc_, pmax_, dmax_, gain_] Produce portrait of encoded information given threshold function of specified gain."

Thrpout::usage="Thrpout[pmin_, pinc_, pmax_, d_, gain_, offset_: 0] Produce portrait of the output probability given threshold functions of specified gain. Offset=0: No tonic inhibition. >0, Floor[d/2] offset."

PDesc::usage="PDesc[f_, TB_] Returns the population acitivty (i.e. probability that a cell in active) implied by frequency 'f' in time bin of TB seconds assuming \
Poisson spike train."

FFromp::usage="FFromp[p_, TB_] Returns the firing rate (in Hz) implied by activity level 'p' assuming that cell is considered active if at least 1 spike occurs \
in timebin of duration TB when firing with Poisson spike train."

(*COMPUTE FUNCTIONS********************************************************************)

Parallelinit::usage = "Parallelinit[]: Distributes definitions in this package across available kernels such that Parallelsubmit[] can be used \
with 'Findent' and 'Entropyapprx'."

FindentNUM::usage = "FindentNUM[tc,\[Mu],d,\[Phi],\[Gamma],pon,reps,external,Ainit,corr]: Numerically determines the entropy of the outputs of a BPN. The channel is used tc times, \
having \[Mu] inputs active with pon, \[Gamma] outputs with d connections each and threshold \[Phi]. The result is averaged over reps network realisations. External=0 causes \
a randomly generated matrix to be used. External=1 allows the argument Ainit to contain the desired matrix. Corr > 0 permits simple correlation structure in MFs \
with spatial extent chosen from exponential distribution with scale parameter 1/corr"

Entvsg::usage = "Entvsg[patts,\[Mu],d,\[Phi],pon,\[Gamma]max,\[Gamma]min,\[Gamma]inc]: Returns an array containing the entropy of the output units of a BPN having \[Mu] inputs \
with d connections per output and output thresholds \[Phi]. patts is the number of events encoded by the inputs, active with probability pon. \
Entropy is found for networks with \[Gamma]min \[Rule] \[Gamma]max outputs in increments of \[Gamma]inc. Makes use of all kernels."

EntvsgDEP::usage="Decpecated version of Entvsg."

EntvsgNUM::usage = "EntvsgNUM[patts,\[Mu],d,\[Phi],pon,\[Gamma]max,\[Gamma]min,\[Gamma]inc,reps]: Numerical version Entvsg, where the result is averaged over reps network realisations. \ 
Returns an array containing the entropy of the output units of a BPN having \[Mu] inputs \
with d connections per output and output thresholds \[Phi]. patts is the number of pattens presented to the inputs, active with probability pon. \
Entropy is found for networks with \[Gamma]min \[Rule] \[Gamma]max outputs in increments of \[Gamma]inc. Makes use of all kernels."

EntSparseSizeNUM::usage = "EntSparseSizeNUM[grd,mfd,\[Mu],pon,d,\[Phi],patts,reps]: Invokes 'FindentNUM' for a range of network sizes and sparsities. Network size is \
fixed with a ratio between number of inputs \[Mu] and number of outputs \[Gamma], which is determined from the granule cell number density, grd and the \
glomerulli density mfd. \[Mu] is a vector of length x specifying the sizes the networks to be considered. pon is a vector of lenght y, of input \
sparsities at which the entropy of these networks should be evaluated. patts is a matrix of dimension x by y, holding the number of randomly \
generated input patterns that are supplied to each network. This function automatically makes use of all available mathematica kernels."

EntSparseSizeDEP::usage = "EntSparseSize[grd,mfd,\[Mu],pon,d,\[Phi],patts]: Invokes 'Findent' for a range of network sizes and sparsities. Network size is \
fixed with a ratio between number of inputs \[Mu] and number of outputs \[Gamma], which is determined from the granule cell number density, grd and the \
glomerulli density mfd. \[Mu] is a vector of length x specifying the sizes the networks to be considered. pon is a vector of lenght y, of input \
sparsities at which the entropy of these networks should be evaluated. patts is a matrix of dimension x by y, holding the number of randomly \
generated input patterns that are supplied to each network. This function automatically makes use of all available mathematica kernels. This function
is deprecated."

Maxeff::usage="Maxeff[fileno, cost, info]: Finds the maximum efficiency and the network attaining that efficiency within the data structure info returned by \ 
EntSparseSize and indexed by the file number, fileno. The cost is the cost for each network as found by Levycost."

InfoPortrait::usage="DoptPortrait[tc_,m_Integer,dini_Integer,dinc_Integer,dmax_Integer,pini_Real,pinc_Real,pmax_Real,g_Integer]: Computes the entropy surface of the outputs of a\
BPN over the number of connections and the probability of input activation. A surface of calculated for each threshold point."
					 
DoptPortraitTFAST::usage="DoptPortrait[tc_,m_Integer,dini_Integer,dinc_Integer,dmax_Integer,pini_Real,pinc_Real,pmax_Real,g_Integer]: Computes the entropy surface of the outputs of a\
BPN over the number of connections and the probability of input activation. A surface of calculated for each threshold point."
					 
DoptPortraitTSLOW::usage="DoptPortrait[tc_,m_Integer,dini_Integer,dinc_Integer,dmax_Integer,pini_Real,pinc_Real,pmax_Real,g_Integer]: Computes the entropy surface of the outputs of a\
BPN over the number of connections and the probability of input activation. A surface of calculated for each threshold point."	

(*UTILITY********************************************************************)

Progress::usage = "Progress[]: Displays a progress bar for EntSparseSize. It is recommended that time consuming jobs are performed on the compute nodes."

Loadfiles::usage="Loadfiles[maxd, fileid, path, fpre, fsuff]: Loads files saved with the Mathematica function Export. maxd is the number of values of d covered by the saved data\
set, fileid is an array of file numbers, fpre is the prefix to those numbers and fsuff is the suffix to those numbers."

LoadPortraits::usage="LoadPortraits[filebase,ext,mlist,pattlist,maxd,pmin,pinc,pmax]: Load data generated by InfoPortrait"

Loaddump::usage="Loaddump[file]: Loads a mathematica variable that has been directly dumped to file with Export[file, <variable_name>]"

Numformat::usage="Numformat[number_, len_]: Returns a string representation of number with 'len' digits to the right of the decimal point. Useful for loading data."

Info::usage = "Info[]: Gives information about this package, including version number."
(**************************************************************************************)

Begin["`Private`"]

(*INTERNAL FUNCTIONS*******************************************************************)
(*surprisal of an event that occurs with probability p*********************************)
h[p_] := If[p == 0, 0, -(Log[2, p])]
(*mean surprisal of an event that occurs with probability p****************************)
h2[p_] := If[p == 0, 0, -p*(Log[2, p])]
(*compute number of outputs given the input (mfdens) and output (grcdens) densities****)
numoutputs[\[Mu]_Integer,mfdens_,grcdens_]:=N[(\[Mu]/mfdens)*grcdens]
(*fast computation of the binomial distribution****************************************)
convbino[N_Integer, p_] := (PDF[BinomialDistribution[N, p], #1] & ) /@ (Range[N + 1] - 1)
r[mu_, sx_] := (Binomial[mu, sx] - 1)/Binomial[mu, sx]
(*fast computation of core calculation within the entropy approxiamtion****************)
Rp = Compile[{
		      {p}, 
	          {sy}, 
	          {pO}, 
              {g},               
              {tjdist},
              {indist}
             }, 
              
              ((If[# <= 0, 1, #]&) 
                              @ (indist*(Plus @@ 
                              Table[((-1)^z*Binomial[g - sy+1,z]*
                              	   (If[(#1 - pO) < 0, 1, If[#1 > 0, (#1 - pO)/#1,1]]
                              	   *(#1^(sy-1 + z)))), 
                              {z, 0, g - sy+1}] + pO &)@p) 
                              + (tjdist/Binomial[g, sy-1]))
            ]
(*optimised version of this function for internal use***********************************)
entropyapprxINTERNAL[tc_,m_Integer,d_Integer,thr_Integer,pin_Real,g_Integer,p_] := 
 Module[{dist = Joint\[Alpha]\[Beta][m, pin, g, d, thr], 
 	     outdist, 
         indist, 
         nc,  
         pO, 
         out},
         
         outdist = Total[Transpose[dist]];
         indist = Total[dist];
         nc = Uniquecodes[indist, m, tc]; 
         pO = Table[If[nc[[i]] > 0, 1/nc[[i]], 0], {i, 1, m + 1}]; 
         (*p  = Table[p\[Phi][m, x, thr, d], {x, 0, m}];*)
          
         out = Plus@@Plus@@Outer[
         	                     dist[[#2]][[#1]]*-1*Log2[Rp[p[[#1]],#2, pO[[#1]], g, Plus@@dist[[#2]][[Complement[Range[m + 1], {#1}]]], indist[[#1]]]] &
         	                    , Range[m+1],Range[g]+1]
         	               -outdist[[1]]*Log2[outdist[[1]]]; 
 out] 
 (*core calculation for determining pattern separation**********************************)
 dk = Compile[{
      		  {inputs},
      		  {pon},
      		  {s1}
      		 },
   				If[pon==0 || pon==1,0,Total[Table[
   												 Total[Table[(s1 + s2 - 2*z)*Binomial[inputs - s1, s2 - z]*
             												 Binomial[s1, z]*pon^s2*(1 - pon)^(inputs - s2), 
           													{z, 0, Min[s1, s2]}]], 
           										 {s2, 0, inputs}]
           								   ]
           		  ]
      		]           
      		
(*PACKAGE FUNCTIONS********************************************************************)
Bitent[p_] := If[p == 0, 0, If[p == 1, 0, -(p*Log[2, p] + (1 - p)*Log[2, 1 - p])]]
(**************************************************************************************)
(*NB This looks subtley different to the C implementation, specifically in the rounding of the unique codes*)
Sourceent[\[Mu]_Integer,pin_,uc_]:=If[pin == 0 || pin == 1,
														 0,
	                                                     Total[
	                                                     	   uc*Map[
	                                                     	   	      h2,Table[If[uc[[s+1]]>=1,
	                                                     	   	      	                       (Binomial[\[Mu],s]/uc[[s+1]]),
	                                                     	   	      	                       If[Round[uc[[s+1]]]>0,
	                                                     	   	      	                       	                     (Binomial[\[Mu],s]/Round[uc[[s+1]]]),0
	                                                     	   	      	                       	       ]
	                                                     	   	      	          ]*(pin^s)*((1-pin)^(\[Mu]-s)),
	                                                     	   	      	       {s,0,\[Mu]}]]]]
(**************************************************************************************)
Pout[d_, pin_, \[Phi]_] := Total[Table[PDF[BinomialDistribution[d, pin],b], {b, \[Phi], d}]]
(**************************************************************************************)
p\[Phi][\[Mu]_Integer, \[Beta]_Integer, \[Phi]_Integer, d_Integer] := Fold[Plus, 0, (PDF[HypergeometricDistribution[d, \[Beta], \[Mu]], #1] & ) /@(Range[If[\[Beta] < d, \[Beta], d] - (\[Phi] - 1)] + (\[Phi] - 1))] 
 (*************************************************************************************)
Typset[r_, H_, n_] := (r/H)*(2^(n*H))
(**************************************************************************************)
Joint\[Alpha]\[Beta][m_Integer, pon_Real, g_Integer,d_Integer, thr_Integer] := 
 Outer[#1*#2 &, Table[1, {g + 1}], (PDF[BinomialDistribution[m, pon], #1] &) /@ (Range[m + 1] - 1)]*
 Outer[(PDF[BinomialDistribution[g, #2], #1] &), (Range[g + 1] - 1), (p\[Phi][m, #1, thr, d] &) /@ (Range[m + 1] - 1)]
(**************************************************************************************)
Uniquecodes[indist_, \[Mu]_Integer, nc_] := Module[{bin = (Binomial[\[Mu], #1] & ) /@ (Range[\[Mu] + 1] - 1)}, 
	                                              N[bin*(1 - Exp[-(Round[indist*nc])/bin]),Dimensions[IntegerDigits[2^\[Mu]]][[1]]]]
(**************************************************************************************)
connectedMscheme[\[Mu]_Integer, d_Integer, draw_Integer, \[Epsilon]_, mdims_Integer] := 
  Table[
  	    Table[If[i == j, 
  	    	            If[i != 1, 
  	    	            	     If[i == mdims, 
  	    	            	     	          1 - \[Epsilon], (d + (i - 1) - (draw - 1))/(\[Mu] - (draw - 1)) - \[Epsilon]],
  	    	            	     	          (d + (i - 1) - (draw - 1))/(\[Mu] - (draw - 1))
  	    	            	        ], 
                                 If[j == i + 1, 
                                 	           (\[Mu] - (d + (i - 1)))/(\[Mu] - (draw - 1)), 
                                 	           If[i == j + 1, 
                                 	           	             \[Epsilon], 
                                 	           	             0
                                 	           	  ]
                                   ]
                ], {j, 1, mdims}], 
  {i, 1, mdims}]
(**************************************************************************************)
Occupancies[nmat_, ini_, \[Gamma]_Integer] := Prepend[(ini . MatrixPower[nmat, #1] & ) /@ Range[\[Gamma] - 1], ini]
(**************************************************************************************)
(*NB: This function is deprecated
A[\[Mu]_, d_, \[Phi]_, \[Gamma]_, sy_, p\[CapitalOmega]_, nc_, indist_] := 
  Table[Table[(If[#1 < 0, 0, #1] & )[
     If[x == y, (Plus @@ Table[(-1)^z*Binomial[\[Gamma]-sy,z]*(If[(#1-p\[CapitalOmega][[x]])<0,1,If[#1>0,(#1-p\[CapitalOmega][[x]])/#1,1]]*(#1^(sy+z))),
             {z,0,\[Gamma]-sy}]+p\[CapitalOmega][[x]] & )@  p\[Phi][\[Mu], x - 1, \[Phi], d], (If[#1 == 1 && sy==\[Gamma],1,(#1^(sy))*((1 - #1)^(\[Gamma]-sy))] & )[
       p\[Phi][\[Mu], x - 1, \[Phi], d]]]], {y, 2, \[Mu] + 1}], {x, 2, \[Mu] + 1}]*)

(*NB: This function is deprecated
EntropyapprxDEP[tc_,\[Mu]_Integer,d_Integer,\[Phi]_Integer,pon_Real,g_Integer]:=Module[{distc,outdist,indist,nc,p\[CapitalOmega],out,y},distc=Joint\[Alpha]\[Beta][\[Mu],pon,g,d,\[Phi]];indist=Total[distc];outdist=Total[Transpose[distc]];
             nc=Uniquecodes[indist,\[Mu],tc];p\[CapitalOmega]=Table[If[nc[[b]]>0,1/nc[[b]],0],{b,1,\[Mu]+1}];
               out=(outdist[[1]]*h[outdist[[1]]]);For[y=2,y<=g+1,y++,If[outdist[[y]]>=0.0001,out=out+(distc[[y]][[2;;\[Mu]+1]]).(h[#]& /@ (indist[[2;;\[Mu]+1]].A[\[Mu],d,\[Phi],g,y-1,p\[CapitalOmega],nc,indist])),out=out]];out]*)
               
(*NB: This function is deprecated
EntvsgDEP[patts_, \[Mu]_Integer, d_Integer, \[Phi]_Integer, pon_, \[Gamma]max_Integer, \[Gamma]min_Integer, \[Gamma]inc_Integer] := 
  Module[{}, DistributeDefinitions[EntropyapprxDEP,A]; Parallelinit[]; WaitAll[Table[ParallelSubmit[{\[Gamma]},EntropyapprxDEP[patts, \[Mu], d, \[Phi],pon,\[Gamma]]], {\[Gamma], \[Gamma]min, \[Gamma]max, \[Gamma]inc}]]]*)               
(**************************************************************************************)
EntropyApprx[tc_,m_Integer,d_Integer,thr_Integer,pin_Real,g_Integer] := 
 Module[{dist = Joint\[Alpha]\[Beta][m, pin, g, d, thr], 
 	     outdist, 
         indist, 
         nc, 
         p, 
         pO, 
         out},
         
         outdist = Total[Transpose[dist]];
         indist = Total[dist];
         nc = Uniquecodes[indist, m, tc]; 
         pO = Table[If[nc[[i]] > 0, 1/nc[[i]], 0], {i, 1, m + 1}]; 
         p  = Table[p\[Phi][m, x, thr, d], {x, 0, m}];
          
         out = Plus@@Plus@@Outer[
         	                     dist[[#2]][[#1]]*-1*Log2[Rp[p[[#1]],#2, pO[[#1]], g, Plus@@dist[[#2]][[Complement[Range[m + 1], {#1}]]], indist[[#1]]]] &
         	                    , Range[m+1],Range[g]+1]
         	               -outdist[[1]]*Log2[outdist[[1]]]; 
 out]    
 (**************************************************************************************)   
(*Note that EntSparseSize and Entvsg are not fully optimised. The table of Entropyapprx calls recalculates the joint number distribution, many times. Depricated
EntSparseSizeDEP[grd_, mfd_, \[Mu]_, pin_, d_Integer, \[Phi]_Integer,patts_] := Module[{m,g, \[Mu]dim = Dimensions[\[Mu]][[1]], pindim = Dimensions[pin][[1]],out}, 
  prog = 0;out = Table[Table[0, {x, 1, pindim}], {y, 1, \[Mu]dim}];For[m = 1, m <= \[Mu]dim, m++,g = Round[numoutputs[\[Mu][[m]], mfd, grd]];g=If[g<1,1,g]; 
         DistributeDefinitions[EntropyapprxDEP, Uniquecodes, Joint\[Alpha]\[Beta], convbino,p\[Phi],A,h]; 
            out[[m]] = WaitAll[Table[ParallelSubmit[{s, m, g},EntropyapprxDEP[patts[[m]][[s]], \[Mu][[m]], d, \[Phi], pin[[s]], g]], {s, pindim}]]; 
            prog = If[\[Mu]dim > 1, (m - 1)/(\[Mu]dim - 1), m - 1/1]];out]*)
(**************************************************************************************)
InfoPortrait[tc_,m_Integer,dini_Integer,dinc_Integer,dmax_Integer,pini_Real,pinc_Real,pmax_Real,g_Integer] := 
 Module[{entropy,
 		dlist,
 		tlist,
 		plist,
 		pvals,
 		listl},
 		
 		 tlist = Flatten[Table[Table[Range[d],{d,dini,dmax,dinc}],{z,pini,pmax,pinc}]];
         dlist = Flatten[Table[Table[Table[d,{z1,d}],{d,dini,dmax,dinc}],{z2,pini,pmax,pinc}]];
         plist = Flatten[Table[Table[Table[z2,{z1,d}],{d,dini,dmax,dinc}],{z2,pini,pmax,pinc}]];
         listl = Dimensions[plist][[1]];
         pvals = Table[Table[p\[Phi][m, x, tlist[[dex]], dlist[[dex]]], {x, 0, m}],{dex,listl}];
 		
         Parallelinit[];
         entropy = WaitAll[
         	         Table[ParallelSubmit[{dex,tlist,plist,dlist,pvals},
         	           	                    entropyapprxINTERNAL[tc,m,dlist[[dex]],tlist[[dex]],plist[[dex]],g,pvals[[dex]]]
         	           	                 ],
         	         {dex,listl}]];
                   entropy
 ]
 (**************************************************************************************)
 LoadPortraits[filebase_,ext_,mlist_,pattlist_,maxd_,pmin_,pinc_,pmax_]:=
 Module[{rawdata,
         outdata,
         pdex,
         mdex,
         pdims=Dimensions[pattlist][[1]],
         mdims=Dimensions[mlist][[1]],
         ppoints,
         netpoints,
         md},
         
         outdata=Table[Table[0,{j,pdims}],{i,mdims}];
         ppoints=Floor[(pmax-pmin)/pinc]+1;
         
         For[mdex=1,mdex<=mdims,mdex++,
         	If[mlist[[mdex]]<maxd,md=mlist[[mdex]],md=maxd];
         	netpoints=Plus@@Range[md];
            For[pdex=1,pdex<=pdims,pdex++,
               rawdata=
                      Flatten[Import[filebase<>ToString[pattlist[[pdex]]]<>"_"<>ToString[mlist[[mdex]]]<>ext]];
               outdata[[mdex]][[pdex]]=
                      Table[Table[Table[rawdata[[(pdex-1)*netpoints+Total[Range[d-1]]+thr]],{pdex,ppoints}],{thr,d}],{d,1,md}];      
            ];
         ];	
 outdata													              
]
(**************************************************************************************) 
(*NB Nmerical functions were not fully optimised*)
EntSparseSizeNUM[grd_, mfd_, \[Mu]_, pin_, d_Integer, \[Phi]_Integer, patts_, reps_] := 
                                      Module[{m, s, g, x, \[Mu]dim = Dimensions[\[Mu]][[1]], pindim = Dimensions[pin][[1]], out}, 
                                              prog = 0; out = Table[Table[0, {x, 1, pindim}], {y, 1, \[Mu]dim}]; 
                                              For[m = 1, m <= \[Mu]dim, m++, g = Round[numoutputs[\[Mu][[m]], mfd, grd]]; 
                                              Parallelinit[]; (*DistributeDefinitions[Binmat, FindentNUM];*) 
                                              out[[m]] = WaitAll[Table[ParallelSubmit[{s, m, g}, FindentNUM[patts[[m]][[s]], \[Mu][[m]], d, \[Phi], g, pin[[s]], reps]], {s, pindim}]]; 
                                              prog = If[\[Mu]dim > 1, (m - 1)/(\[Mu]dim - 1), m - 1/1]]; 
                                             out]
(**************************************************************************************)
Parallelinit[] := DistributeDefinitions[entropyapprxINTERNAL,
										Entropyapprx, 
										Uniquecodes, 
										Joint\[Alpha]\[Beta], 
										p\[Phi],
										h,
										FindentNUM,
										Binmat,
										Rp,
										Hdist,
										dk,
										Pout]
(**************************************************************************************)
Progress[]:=ProgressIndicator[Dynamic[prog]]
(**************************************************************************************)
Entvsg[patts_, \[Mu]_Integer, d_Integer, \[Phi]_Integer, pon_, \[Gamma]max_Integer, \[Gamma]min_Integer, \[Gamma]inc_Integer] := 
  Module[{}, 
  	     Parallelinit[]; 
  	     WaitAll[Table[
  	     	           ParallelSubmit[{\[Gamma]},EntropyApprx[patts, \[Mu], d, \[Phi],pon,\[Gamma]]], 
  	     	          {\[Gamma], \[Gamma]min, \[Gamma]max, \[Gamma]inc}]
  	     	    ]
  	     ]
(**************************************************************************************)
Binmat[done_, dtwo_, s_] := 
Module[{rs, n, out = Table[Table[0, {m, 1, dtwo}], {n, 1, done}]}, 
   For[n = 1, n <= done, n++, rs = RandomSample[Range[dtwo], s]; 
                              out[[n]] = Table[If[MemberQ[rs, x], 1, 0], {x, 1, dtwo}]]; 
   out]
 (**************************************************************************************)     
Hdist[\[Mu]_, dmax_, \[Gamma]_, pini_, pinc_, pmax_, mode_] := 
 Module[{distance,	
 		 ptab,
 		 pdftab,
 		 tablen2=Total[Range[dmax]]*(Round[(pmax-pini)/pinc]+1)*(\[Gamma]+1),
 		 dlist,
 		 slist},
   	
 	Switch[mode,1,
 	 
     distance = Total[Table[
     	                   WaitAll[Table[
     	                   	     ParallelSubmit[{pon,s1},
     	                   	     	           PDF[BinomialDistribution[\[Mu], pon], s1]*dk[\[Mu], pon, s1]],
     	                   {pon,pini,pmax,pinc}]],
     	              {s1,0,\[Mu]}]
     	              ]
                
               ,2, 
     
     Parallelinit[];
     slist	  = Table[1,{s2,0,\[Gamma]}];
     ptab	  = Flatten[Table[Table[WaitAll[Table[ParallelSubmit[{d,pon,th},slist*Pout[d,pon,th]],{pon,pini,pmax,pinc}]],{th,d}],{d,dmax}]];
     pdftab   = Flatten[Table[
     	             Table[
     	             	  WaitAll[Table[ParallelSubmit[{pon,d,th},
     	             	  	   Total[
     	             	  	   	    Transpose[
     	             	  	   	        Joint\[Alpha]\[Beta][\[Mu],pon,\[Gamma],d,th]
     	             	  	   	             ]
     	             	  	   	    ]],
     	             	  {pon,pini,pmax,pinc}]],
     	             {th,d}],
     	        {d,dmax}]];
     	        
     dlist    = WaitAll[Table[ParallelSubmit[{ptab,pdftab,dex},(pdftab[[dex]]*dk[\[Gamma],ptab[[dex]],Mod[dex-1,\[Gamma]+1]])], {dex,tablen2}]];
     dlist    = Total[Transpose[Partition[dlist,\[Gamma]+1]]];  
     distance = Table[Table[Table[dlist[[((Total[Range[d-1]]+th)-1)*(Round[(pmax-pini)/pinc]+1)+Round[(pon-pini)/pinc]+1]],{pon,pini,pmax,pinc}],{th,d}],{d,dmax}]
     	                      
               ,3,
              
     Parallelinit[];
     slist	  = Table[1,{s2,0,\[Gamma]}];
     ptab	  = Flatten[Table[Table[WaitAll[Table[ParallelSubmit[{d,pon,th},slist*Pout[d,pon,th]],{pon,pini,pmax,pinc}]],{th,d}],{d,dmax}]];
     pdftab   = WaitAll[Table[ParallelSubmit[{ptab,dex},PDF[BinomialDistribution[\[Gamma], ptab[[dex]]], Mod[dex-1,\[Gamma]+1]]],{dex,tablen2}]];     	        
     dlist    = WaitAll[Table[ParallelSubmit[{ptab,pdftab,dex},(pdftab[[dex]]*dk[\[Gamma],ptab[[dex]],Mod[dex-1,\[Gamma]+1]])], {dex,tablen2}]];
     dlist    = Total[Transpose[Partition[dlist,\[Gamma]+1]]];  
     distance = Table[Table[Table[dlist[[((Total[Range[d-1]]+th)-1)*(Round[(pmax-pini)/pinc]+1)+Round[(pon-pini)/pinc]+1]],{pon,pini,pmax,pinc}],{th,d}],{d,dmax}]
                                         
 	];
 	          
distance]
(**************************************************************************************)
HdistNUM[\[Mu]_, d_, \[Phi]_, \[Gamma]_, pon_, pairs_] := 
   Module[{p, A, in1, in2, out1, out2, indist = 0, outdist = 0}, 
         For[p = 1, p <= pairs, p++, 
            A = Transpose[Binmat[\[Gamma], \[Mu], d]]; 
            in1 = Table[If[RandomReal[] <= pon, 1, 0], {x, 1, \[Mu]}]; 
            in2 = Table[If[RandomReal[] <= pon, 1, 0], {x, 1, \[Mu]}]; 
            indist = indist + Total[Abs[in1 - in2]]; 
            out1 = Thread[Boole[Thread[in1 . A >= \[Phi]]]]; 
            out2 = Thread[Boole[Thread[in2 . A >= \[Phi]]]]; 
            outdist = outdist + Total[Abs[out1 - out2]]
            ]; 
         {outdist,indist}]
(**************************************************************************************)
FindentNUM[tc_, \[Mu]_, d_, \[Phi]_, \[Gamma]_, pon_, reps_, external_: 0, Ainit_: 0, corr_: 0] := 
  Module[{r, A, ent = 0, tmp, rvals, inputs, rnd}, 
                                            For[r = 1, r <= reps, r++, 
                                                                A = Switch[external, 
                                                                	               0, 
                                                                	               Transpose[Binmat[\[Gamma], \[Mu], d]], 
                                                                	               1, 
                                                                	               Ainit];
                                                                tmp = Round[RandomInteger[PoissonDistribution[If[corr==0,1,1/corr]], \[Mu]] + 1];
                                                                rvals = tmp[[1 ;; Min[Flatten[Position[Thread[Accumulate[tmp] >= \[Mu]], True]]]]];
                                                                inputs = If[corr <= 0, 
                                                                	                Table[Table[If[RandomReal[] <= pon, 1, 0], {x, 1, \[Mu]}], {x, 1,tc}], 
                                                                	                Table[Flatten[Map[(Table[rnd, {x, #}] /. rnd -> If[RandomReal[] <= pon, 1, 0]) &, rvals]][[1 ;; \[Mu]]], {x, 1, tc}]];
                                                                ent = ent + N[Entropy[2, Boole[(Thread[#1 >= \[Phi]] &) /@ 
         Table[inputs[[x]].A, {x, 1, tc}]]]]];
  ent/reps]
 (**************************************************************************************)            
EntvsgNUM[patts_, \[Mu]_Integer, d_Integer, \[Phi]_Integer, pon_, \[Gamma]max_Integer, \[Gamma]min_Integer, \[Gamma]inc_Integer,reps_Integer] := 
  Module[{}, 
  	     Parallelinit[]; 
  	     WaitAll[Table[ParallelSubmit[{\[Gamma]}, 
  	     	                          FindentNUM[patts, \[Mu], d, \[Phi], \[Gamma], pon,reps]], 
  	     	          {\[Gamma], \[Gamma]min, \[Gamma]max, \[Gamma]inc}]]]
(**************************************************************************************)
Levycost[\[Mu]_Integer, pin_, x_, \[Gamma]_Integer, pout_, y_,\[Delta]_] := 
  \[Mu]*(pin*x + (1 - pin)) + \[Delta]*\[Gamma]*(pout*y + (1 - pout))
(**************************************************************************************)
Levycost2[\[Mu]_Integer, pin_, x_, \[Gamma]_Integer, pout_, y_,\[Delta]_,G_,d_,connscal_] := 
  (\[Mu]*(pin*x + (1 - pin))) + (\[Delta]*\[Gamma]*(pout*y + (1 - pout)))+(\[Gamma]*d*pin*G*connscal/\[Mu])
(**************************************************************************************)
Maxeff[fileno_, cost_, info_] := Module[{meff, dims1 = Dimensions[info], eff}, 
   eff = Table[Table[Flatten[info[[fileno]][[d]][[\[Phi]]], 1]/cost[[fileno]][[d]][[\[Phi]]], {\[Phi], 1, d}], {d, 1, dims1[[2]]}]; meff = Max[eff]; {meff, Position[eff, meff]}]
(**************************************************************************************)
Loadfiles[maxd_, fileid_, path_, fpre_, fsuff_] := 
   Module[{out,f,file,import,d},
   	      If[fileid=={},out=Table[0, {y, maxd}];
   	      file = StringJoin[path, fpre, fsuff];
   	      import=ReadList[file, "String"];
   	      For[d = 1, d <= maxd, d++, out[[d]] = ToExpression[StringSplit[import[[d]], "\t"]]],
   	                                 out=Table[Table[0, {y, 1, maxd}], {x, 1, Dimensions[fileid][[1]]}]; 
   	                                 For[f = 1, f <= Dimensions[fileid][[1]], f++, 
   	                                 	 file = StringJoin[path, fpre, ToString[fileid[[f]]], fsuff]; 
   	                                 	 import = ReadList[file, "String"]; 
                                         For[d = 1, d <= maxd, d++, out[[f]][[d]] = ToExpression[StringSplit[import[[d]], "\t"]]]]]; 
                                         out]
(**************************************************************************************)
GThreshold[pin_, d_, gain_, offset_: 0] := 
Module[{intmems}, 
  		intmems = IntervalMemberQ[Table[Interval[{(x - 1)*1/(gain*d), If[x == d, 1, x*1/(gain*d)]}], {x, d}], pin]; 
        Thread[If[Thread[# > d], d, #]] &@(Last[Position[intmems, True]][[1]] + If[offset>0,Floor[(d+1)/2],1]-1)];
(**************************************************************************************)
Thrent[data_, pmin_, pinc_, pmax_, d_, gain_, offset_: 0] := 
Table[
      data[[1]][[d]][[GThreshold[pin, d, gain, offset]]][[1]][[1]][[1 + Floor[(pin - pmin)/pinc]]], 
     {pin, pmin, pmax, pinc}]
 (**************************************************************************************)
Thrprune[data_,pmin_, pinc_, pmax_, dmax_, gain_] :=
   {Table[
        Table[
        	 {{Table[
                    If[#!=\[Phi],0,data[[1]][[d]][[#]][[1]][[1]][[1 + Floor[(pin - pmin)/pinc]]]]&@GThreshold[pin, d, gain],
         	 {pin, pmin, pmax, pinc}]}}, 
        {\[Phi], d}],
   {d, dmax}]}
(**************************************************************************************)
Thrpout[pmin_, pinc_, pmax_, d_, gain_, offset_: 0] := Table[Pout[d, pin, GThreshold[pin, d, gain, offset]], {pin, pmin, pmax, pinc}]                                         
(**************************************************************************************)        
Loaddump[file_]:=
                      Flatten[
			 			ToExpression[ 
		      			 StringSplit[
			   			  ReadList[
				        		   file, 
				                   "String"], 
			   		     "\t"]], 
            		   1] 
(**************************************************************************************)
Numformat[number_, len_] := 
              Module[{lhs = RealDigits[number][[2]], 
                      rhs = StringLength[ToString[number - Floor[number]]] - 2}, 
                      If[lhs < 1, lhs = 0, Null];
                      StringDrop[ToString[PaddedForm[number, {lhs + rhs, len}]],1]];     
(**************************************************************************************)   
PDesc[f_, TB_:(3*10^-2)] := N[1 - PDF[PoissonDistribution[f*TB], 0]]
(**************************************************************************************)
FFromp[p_, TB_:(3*10^-2)] := -Log[1 - p]/TB
(**************************************************************************************)
Info[] := Print["BinaryPartitionNetworks is a Mathematica package for working with Binary Partition Networks. Guy \
Billings UCL 2010. v0.5."]                                       
(**************************************************************************************)
End[]
Info[]
EndPackage[]

