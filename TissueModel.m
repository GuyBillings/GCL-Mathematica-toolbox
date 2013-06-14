(* Mathematica Package *)

(* Created by the Wolfram Workbench Feb 14, 2012 *)
(* This package contains functions for setting up and investigating the Anatomical GCL model *)

BeginPackage["TissueModel`"]
(* Exported symbols added here with SymbolName::usage *)

grcinvol::usage="grcinvol[dim_, geom_, dlen_: 15] Granule cells volume of type geom (sphere, cube or torus) of internal dimension dim."

glominvol::usage="glominvol[dim_, geom_] Glomeruli in volume type geom (sphere, cube or torus) of internal dimension dim."

glomdens::usage="glomdens[bigdim_, smalldim_, glompos_, geom_] Density of glomeruli in subvolume of type geom of dimensions smalldim"

grcpos::usage="grcpos[bigdim_, smalldim_, cells_, geom_] Randomly place granule cells in a subvolme of type geom of dimensions smalldim.\
smalldim is the diameter if the geom is a sphere."

glopos::usage="glopos[dim_, nfibres_, rawdist_, distmeanX_, distmeanY_, distmeanZ_, geom_] Randomly place glomeruli in volume of type geom,\
such that they are realistically sampled from mosssy fibres"

glomcount::usage="glomcount[bigdim_, smalldim_, glompos_, geom_] Count the glomeruli in the specified volume given the positions, glompos."

confineglom::usage="confineglom[gloms_, geom_, slen_, dlen_: 15] Wrapper that confines the positions of the glomeruli to a given subvolume"

ALGconnections::usage="ALGconnections[dlen_, mflist_, grclist_, d_, mossyrepeat_, force_] Algorithmic connections function : 
Ensures that the connections have an mf degree distribution that \
matches the random connection matrix and that all granule cells have exactly d connections.\ 
Aims to give an average connection length that is around a target value, \ 
but this is not gauranteed and the resultant distribution will \
depend on other parameters and sampling geometry. \
mossyrepeat = 0 : prevents TARGET degree distribution from reconnecting to identical mossy fibres through independent connections,\
mossyrepeat = 1 : Allows repeat sampling of mossy fibres.\
mossyrepeat = 2 : does not allow repeat sampling and also forces the mossy fibres to be connected with equal probability\
(i.e. does not randomly choose from glom list.)\
mossyrepeat also changes the degree distribution process so that the degree distribution is for granule cells that can only \
connect to a single mossy fibre once.\
force = 0 : Does not force a match between DDist and the degree distribution in the connection list.\
force = 1 : Force the degree distribution to be matched by the model connectivity, by rearranging connections."

grclist2DDist::usage="grclist2DDist[grclist_] Convert the listing of granule cells connected to each glomeruli into a degree distribution."

grclist2conns::usage="grclist2conns[grclist_, glompos_] Convers the listing of granule cells connected to each glomeruli into a connection list."

uniquemfsFROMgrclist::usage="uniquemfsFROMgrclist[grclist_] Degree distribution counting unique connections between granules and mossies."

grcrepeats::usage="grcrepeats[conns_] Counts the number of granule cells that do not connect to d distinct mossy fibres"

netmossycount::usage="netmossycount[conns_] Number of mossy fibres connected to the network"

denlendist::usage="denlendist[conns_, glompos_, granpos_] Distribution of granule cell dendrite length."

meandens::usage="Average numerical density of glomeruli given realisation of some number of mossy fibres in the space."

targetmf::usage="targetmf[targetglom_, biglen_, smalllen_, rdist_, glox_, gloy_, gloz_,geom1_, geom2_, reps_] Determine number of mossy fibres for target density"

corrvsdist::usage="corrvsdist[shellinc_, shellmax_, nets_, rdist_, slen_, centrallen_, Nfib_, Ngrc_, glox_, gloy_, gloz_, denlen_, d_, pin_, \[Phi]_, tc_, geom1_, geom2_, mode_, d\[Mu]_: 15, d\[Sigma]_: 5] \
Covariance between granule cells as one travels radial across a volume of tissue"
 
meangloms::usage="meangloms[nets_, rdist_, slen_, clen_, Nfib_, Ngrc_, glox_, gloy_, gloz_, denlen_, d_, geom1_, geom2_, mode_, d\[Mu]_: 15, d\[Sigma]_: 5, repeat_, force_] Mean number of glomeruli in tissue volume"

Amatrix::usage="Amatrix[conns_] Adjacency matrix formed from the connection list."

neighbourhood::usage="neighbourhood[A_] Computes the neighbourhood matrix of adjacency matrix A"

degreedist::usage="degreedist[nets_, rdist_, slen_, clen_, Nfib_, Ngrc_, glox_, gloy_, gloz_, denlen_, d_, mode_, geom1_, geom2_, algo_: 1, d\[Mu]_: 15, d\[Sigma]_: 5, mossyrepeat_, force_] \
Determine the degree distibution of the glomeruli"
                   
anatomicalent::usage="anatomicalent[nets_, slen_, Nfib_, rdist_, glox_, gloy_, gloz_, geom1_, centrallen_, Ngrc_, geom2_, denlen_, d_, mode_, d\[Mu]_, d\[Sigma]_, tc_, \[Phi]_, pin_, reps_, algo_, repeats_, force_]\ 
Entropy of binary model approximating the anatomical model"
   
tissueprof::usage="tissueprof[rdist_, bigdim_, cendim_, nmfs_, dlen_, d_, glox_, gloy_, gloz_, geom1_, geom2_, mode_, force_] Function to scan these variables across network parameters.\ 
output : {{number mf}, {number mf connected}, {number grc reconnected to same mf}, {number unshared}, {number connected glomeruli},\
{number of granule cells with fewer than target connections}, {number grcs connecting same glomeruli more than once},{distribution of unique connections per cell}}"

meanconpro::usage="meanconpro[rdist_, bigdim_, cendim_, nmfs_, dlen_, d_, glox_, gloy_, gloz_, geom1_, geom2_, reps_, mode_, force_] \
Function to determine the average number of cells that fail to make d unique connections and the average number of connections to unique mossy fibres per cell"                   
                   
Begin["`Private`"]
(* Implementation of the package *)
(*INTERNAL FUNCTIONS*******************************************************************)
glomonfibre[cdf_] := (Position[#, Min[#]] &@ Abs[RandomReal[] - cdf])[[1]][[1]] + 1
memberContact[origin_, point_, dlen_] := If[# <= dlen, #, 0] &@(Sqrt[(#1^2 + #2^2 + #3^2)] & @@ (point - origin)) 
(*Determine if point is within specified Torus, at origin with major radius in the xy plane*) 
torusContact[majorR_, minorr_, point_] := If[# <= minorr, 1, 0] &@(Sqrt[(majorR - Sqrt[#1^2 + #2^2])^2 + #3^2]) & @@ point
(*Return the mossy fibre and glomerulus indicies of the glomerulus that is nearest to the target distance*)
closestallowedMF[attachedmfs_, grccoords_, mflist_, dlen_] := Module[{ad, \[Mu] = Dimensions[mflist][[1]]}, 
                                                                      ad = Table[
                                                                      	         Table[
                                                                      	         	    If[MemberQ[attachedmfs, {m, l}], \[Infinity], Abs[dlen - EuclideanDistance[grccoords, mflist[[m]][[l]]]]], 
                                                                      	         	   {l,Dimensions[mflist[[m]]][[1]]}], 
                                                                      	         {m, \[Mu]}];
                                                                      Position[ad, Min[ad]][[1]]
                                                                     ]
                                                                     closestallowedGRC[attachedgrcs_, grclist_, mfcoords_, connlist_, d_, dlen_, omit_] := 
                  Module[{allowed1, allowed2, out, ad, full, \[Gamma] = Dimensions[grclist][[1]]}, 
                          allowed1 = Complement[Range[\[Gamma]], attachedgrcs]; 
                          full = DeleteCases[Table[If[# >= d, g, Null] &@Dimensions[connlist[[g]]][[1]], {g, \[Gamma]}], Null]; 
                          allowed2 = Complement[Range[\[Gamma]], full]; 
                          out = Intersection[allowed1, allowed2, Complement[Range[\[Gamma]], omit]]; 
                          ad = Abs[dlen - Thread[EuclideanDistance[mfcoords, #]] & /@ grclist]; 
                          {If[Dimensions[#][[1]] != 0, #[[1]][[1]], 0] &@ Position[ad, Min[ad[[out]]]], out, allowed1, allowed2}
                        ]                        
(*Function to shuffle connections to allow perfect match between specified degree distribution and algorithmic connections*)
shuffleconns[granpos_, glompos_, d_, dlen_, incompletes_, grcs_, DDist_, connlist_, mode_] := 
             Module[{nincompletes = Dimensions[incompletes][[1]], target, allowed,
                     numbertoshuffle, newconnlist = connlist, swap, bestswap, 
                     \[Gamma] = Dimensions[granpos][[1]], grcsL = grcs}, 
                     For[target = 1, target <= nincompletes, target++, 
                        numbertoshuffle = DDist[[incompletes[[target]][[1]]]][[incompletes[[target]][[2]]]] - 
                                          Dimensions[grcs[[incompletes[[target]][[1]]]][[incompletes[[target]][[2]]]]][[1]]; 
                        For[swap = 1, swap <= numbertoshuffle, swap++, 
                            allowed = closestallowedGRC[{}, granpos, glompos[[incompletes[[target]][[1]]]][[incompletes[[target]][[2]]]], newconnlist, d, dlen, {}][[2]]; 
                            bestswap = optswap[granpos, glompos, grcsL, allowed, Switch[mode, 1, 
                                                                                        Complement[Range[\[Gamma]], Union[allowed, Flatten[grcsL[[incompletes[[target]][[1]]]][[incompletes[[target]][[2]]]]]]], 
                                                                                        0, 
                                                                                        Complement[Range[\[Gamma]], Union[allowed, Flatten[grcsL[[incompletes[[target]][[1]]]]]]],
                                                                                        2, 
                                                                                        Complement[Range[\[Gamma]], Union[allowed, Flatten[grcsL[[incompletes[[target]][[1]]]]]]]
                                                                                       ], incompletes[[target]], newconnlist, dlen, mode];
                            newconnlist[[bestswap[[2]]]] = DeleteCases[newconnlist[[bestswap[[2]]]], bestswap[[3]]]; 
                            grcsL[[bestswap[[3]][[1]]]][[bestswap[[3]][[2]]]] = DeleteCases[grcsL[[bestswap[[3]][[1]]]][[bestswap[[3]][[2]]]], bestswap[[2]]];
                            AppendTo[newconnlist[[bestswap[[2]]]], incompletes[[target]]]; 
                            AppendTo[grcsL[[incompletes[[target]][[1]]]][[incompletes[[target]][[2]]]], bestswap[[2]]]; 
                            AppendTo[newconnlist[[bestswap[[1]]]], bestswap[[3]]]; 
                            AppendTo[grcsL[[bestswap[[3]][[1]]]][[bestswap[[3]][[2]]]], bestswap[[1]]]
                           ]
                       ];
                     {newconnlist, grcsL}]
(*Swap GC connections to satisfy connectivity constraints*)                     
(*The swap process: 1)Cell A -> disconnect incomplete glom (i), 2)Cell B -> connect to i, 3)Cell B -> disconnect from another 
complete mf (c) not yet connected to Cell A, 4) Cell A -> connect to c, 5) Cell A -> Reconnect to i. This simplfies to steps 2 through 4. 
Data structure: {Cell A, Cell B,{c:m,c:gl}}. The swap that gives the minimal mean squared deviation (across connections involved) from the 
target is output*) 
(*For mossyrepeat = 0 this function must omit those CellA that are connected already to the complete glomerulus to ensure that step 4 
does not introduce reconnections to the same mf.*)
optswap[granpos_, glompos_, grcs_, CellAcands_, CellBcands_, incomglom_, connlist_, dlen_, mode_] := 
        Module[{rawdist, comglom, optpos, omit = Table[{}, {x, Dimensions[CellBcands][[1]]}]}, 
               comglom = Table[Table[Switch[mode, 1, 
               	                            connlist[[CellBcands[[x]]]], 
               	                            0, 
               	                            Complement[connlist[[CellBcands[[x]]]], connlist[[CellAcands[[y]]]]], 
               	                            2, Complement[connlist[[CellBcands[[x]]]], connlist[[CellAcands[[y]]]]]
               	                           ], 
               	                     {x, Dimensions[CellBcands][[1]]}], 
               	               {y, Dimensions[CellAcands][[1]]}];
               rawdist = Table[
               	               Table[
               	                     Table[
               	                     	   ((EuclideanDistance[granpos[[CellBcands[[cellB]]]], glompos[[incomglom[[1]]]][[incomglom[[2]]]]] - dlen)^2 + 
               	                     	   (EuclideanDistance[granpos[[Complement[CellAcands, omit[[cellB]]][[cellA]]]], glompos[[comglom[[cellA]][[cellB]][[glom]][[1]]]][[comglom[[cellA]][[cellB]][[glom]][[2]]]]] - dlen)^2)/2, 
               	                     	   {glom, Dimensions[comglom[[cellA]][[cellB]]][[1]]}], 
               	                     {cellA,Dimensions[Complement[CellAcands, omit[[cellB]]]][[1]]}], 
               	               {cellB, Dimensions[CellBcands][[1]]}]; 
               optpos = Position[rawdist, Min[rawdist]][[1]]; 
               {CellAcands[[optpos[[2]]]], CellBcands[[optpos[[1]]]], comglom[[optpos[[2]]]][[optpos[[1]]]][[optpos[[3]]]], rawdist, omit}
              ]
(*Compute radial correlation between cells*)               
cellcor[tc_, \[Phi]_, pon_, A_] := Module[{output, \[Mu] = Dimensions[A[[1]]][[2]], cmat},
                                           output = Boole[(Thread[# >= \[Phi]] &) /@ Table[Table[
                                                                                                 If[RandomReal[] <= pon, 1, 0],
                                                                                                 {x, 1, \[Mu]}].Transpose[A[[1]]], 
                                                                                           {x, 1, tc}]
                                                         ]; 
                                           cmat = Covariance[output]; 
                                           cmat
                                         ]
(*Compute distance between cells*)
celldist[granpos_] := Table[Table[EuclideanDistance[granpos[[one]], granpos[[two]]], {one, Dimensions[granpos][[1]]}], {two, Dimensions[granpos][[1]]}]
(*count unique connected GCs*)
uniqueconns[conns_] := Module[{grc, numgrc = Dimensions[conns][[1]], d, cnt = 0}, 
                                For[grc = 1, grc <= numgrc, grc++, 
                                    d = If[# == {}, 0, #[[1]]] &@Dimensions[conns[[grc]]]; 
                                    If[d != 0, cnt = cnt + If[# == 0, 1, 0] &@ Total[Table[Dimensions[Flatten[Position[Flatten[conns[[All, All, 1]][[Complement[Range[numgrc], {grc}]]]], conns[[grc]][[All, 1]][[x]]]]][[1]] - 1, 
                                    	                                                   {x, d}]
                                    	                                            ]
                                      ]
                                   ]; 
                                   cnt
                              ]     
(*PACKAGE FUNCTIONS*******************************************************************)                   
grcinvol[dim_, geom_, dlen_: 15] := Switch[geom, "cube", 
 	                                       Round[1.9*10^6*(dim/10^3)^3], 
 	                                       "sphere", 
 	                                       Round[1.9*10^6*(4/3)*Pi*(dim/2/10^3)^3], 
 	                                       "torus", 
 	                                       Round[1.9*10^6*2*Pi^2*(dim/10^3)*(dlen/10^3)^2]
 	                                      ]
(**************************************************************************************) 	                                
 glominvol[dim_, geom_] := Switch[
 	                              geom, "cube", 
 	                              Round[6.6*10^5*(dim/10^3)^3], 
 	                              "sphere", 
                                  Round[6.6*10^5*(4/3)*Pi*(dim/2/10^3)^3]
                                 ]	                                      
(**************************************************************************************)                
glomdens[bigdim_, smalldim_, glompos_, geom_] := Switch[
	                                                    geom, "cube", 
                                                        N[glomcount[bigdim, smalldim, glompos, geom]/(smalldim/10^3)^3], 
                                                        "sphere", 
                                                        N[glomcount[bigdim, smalldim, glompos, geom]/((4/3)*Pi*(smalldim/2/10^3)^3)]
                                                       ]
(**************************************************************************************)
grcpos[bigdim_, smalldim_, cells_, geom_] := Module[{polar, coords = Table[0, {c, cells}], i}, 
                                                    Switch[geom, "cube", 
                                                           coords = Table[
                                                           	              {RandomReal[{-smalldim/2, smalldim/2}], 
                                                                           RandomReal[{-smalldim/2, smalldim/2}], 
                                                                           RandomReal[{-smalldim/2, smalldim/2}]}, 
                                                                         {x, cells}],
                                                           "sphere", 
                                                           For[i = 1, i <= cells, i++, 
                                                           	   (*WARNING, seems that uniform sampling in spherical coords does not map back to uniform
                                                           	   	sampling in cartesian space*)
                                                               polar = {RandomReal[{0, smalldim/2}], RandomReal[{0, 2*Pi}], RandomReal[{0, 2*Pi}]}; 
                                                               coords[[i]] = {polar[[1]]*Sin[polar[[2]]]*Cos[polar[[3]]], 
                                                               	              polar[[1]]*Sin[polar[[2]]]*Sin[polar[[3]]], 
                                                                              polar[[1]]*Cos[polar[[2]]]}
                                                              ], 
                                                           "torus", 
                                                           For[i = 1, i <= cells, i++, 
                                                               polar = {RandomReal[{0, 2*Pi}], RandomReal[{0, 2*Pi}]}; 
                                                               coords[[i]] = {(bigdim + smalldim*Cos[polar[[2]]])*Cos[polar[[1]]], 
                                                               (bigdim + smalldim*Cos[polar[[2]]])*Sin[polar[[1]]], 
                                                               smalldim*Sin[polar[[2]]]}
                                                              ]
                                                          ]; 
                                                          coords]
(**************************************************************************************)                                                          
glopos[dim_, nfibres_, rawdist_, distmeanX_, distmeanY_, distmeanZ_, geom_] := 
                                                          Module[{primary = Table[{}, {f, nfibres}], x, y, polar, ptemp, pglom, cdfpglom, coords}, 
                                                                 pglom = rawdist/Total[rawdist]; 
                                                                 cdfpglom = Accumulate[pglom]; 
                                                                 Switch[geom, "cube", 
                                                                        primary = Table[Table[If[y == 1, {RandomReal[{-dim/2, dim/2}], RandomReal[{-dim/2, dim/2}], RandomReal[{-dim/2, dim/2}]}, 
                                                                                                         {RandomReal[ExponentialDistribution[1/distmeanX]], RandomReal[ExponentialDistribution[1/distmeanY]], RandomReal[ExponentialDistribution[1/distmeanZ]]}], 
                                                                                              {y, glomonfibre[cdfpglom]}], 
                                                                                        {x, nfibres}]; 
                                                                        coords = Table[Table[primary[[x]][[1]] + 
                                                                        	                 If[y == 1, 0, Total[primary[[x]][[2 ;; y]]]], 
                                                                        	                 {y, 1, Dimensions[primary[[x]]][[1]]}], 
                                                                                       {x, Dimensions[primary][[1]]}], 
                                                                        "sphere", 
                                                                         For[x = 1, x <= nfibres, x++, 
                                                                         	 For[y = 1, y <= glomonfibre[cdfpglom], y++, 
                                                                                 polar = {RandomReal[{0, dim/2}], RandomReal[{0, 2*Pi}], RandomReal[{0, 2*Pi}]}; 
                                                                                 ptemp = If[y == 1, 
                                                                                 	        {polar[[1]]*Sin[polar[[2]]]*Cos[polar[[3]]], polar[[1]]*Sin[polar[[2]]]*Sin[polar[[3]]], polar[[1]]*Cos[polar[[2]]]}, 
                                                                                            {RandomReal[ExponentialDistribution[1/distmeanX]], RandomReal[ExponentialDistribution[1/distmeanY]], RandomReal[ExponentialDistribution[1/distmeanZ]]}
                                                                                           ]; 
                                                                                 AppendTo[primary[[x]], ptemp]]]; 
                                                                                 coords = Table[Table[primary[[x]][[1]] + 
                                                                                 	                  If[y == 1, 0, 
                                                                                 	      	             Total[primary[[x]][[2 ;; y]]]],
                                                                                 	      	          {y, 1, Dimensions[primary[[x]]][[1]]}], 
                                                                                                {x, Dimensions[primary][[1]]}]
                                                                       ]; 
                                                                      coords]; 
(**************************************************************************************)                                                                      
glomcount[bigdim_, smalldim_, glompos_, geom_] := Module[{gloms, cnt}, 
	                                                     gloms = Flatten[glompos, 1]; 
                                                         Switch[geom, "cube", 
                                                                cnt = IntervalMemberQ[Interval[{-smalldim/2, smalldim/2}], #[[1]]] && 
                                                                      IntervalMemberQ[Interval[{-smalldim/2, smalldim/2}], #[[2]]] &&
                                                                      IntervalMemberQ[Interval[{-smalldim/2, smalldim/2}], #[[3]]] & /@ gloms; 
                                                                cnt = Dimensions[Cases[cnt, True]][[1]], 
                                                                "sphere", 
                                                                cnt = IntervalMemberQ[Interval[{0, smalldim/2}], Sqrt[#[[1]]^2 + #[[2]]^2 + #[[3]]^2]] & /@ gloms; 
                                                                cnt = Dimensions[Cases[cnt, True]][[1]]
                                                               ]; 
                                                              cnt]
(**************************************************************************************)                                                              
confineglom[gloms_, geom_, slen_, dlen_: 15] := Module[{out, mf, mfs = Dimensions[gloms][[1]], gls, glo}, 
                                                        out = Table[{}, {x, mfs}]; 
                                                        For[mf = 1, mf <= mfs, mf++, 
                                                        	gls = Dimensions[gloms[[mf]]][[1]]; 
                                                            For[glo = 1, glo <= gls, glo++, 
                                                                Switch[geom, "sphere", 
                                                                       If[# != 0, AppendTo[out[[mf]], gloms[[mf]][[glo]]], Null] &@ memberContact[{0, 0, 0}, gloms[[mf]][[glo]], slen/2], 
                                                                       "cube", 
                                                                       If[# != 0, AppendTo[out[[mf]], gloms[[mf]][[glo]]], Null] &@ If[(Abs[gloms[[mf]][[glo]][[1]]] <= slen/2) && (Abs[gloms[[mf]][[glo]][[2]]] <= slen/2) && (Abs[gloms[[mf]][[glo]][[3]]] <= slen/2), 1, 0], 
                                                                       "torus", 
                                                                       If[# != 0, AppendTo[out[[mf]], gloms[[mf]][[glo]]], Null] &@ torusContact[slen, dlen, gloms[[mf]][[glo]]
                                                                      ]
                                                               ]
                                                            ]
                                                      ]; 
                                                     out]
(**************************************************************************************)
ALGconnections[dlen_, mflist_, grclist_, d_, mossyrepeat_:2, force_:1] := 
               Module[{connection, newstruct, cnt = 0, connlist, glomlist, 
                       DDist, s, staken, g, l, \[Mu] = Dimensions[mflist][[1]], 
                       \[Gamma] = Dimensions[grclist][[1]], connections, m, gloms, 
                       omit = Table[{}, {x, Dimensions[mflist][[1]]}], debug = 0, 
                       gl, grcs, thisgrc, noncomplete}, 
                       connections = \[Gamma]*d; 
                       connlist = Table[{}, {g, \[Gamma]}]; 
                       glomlist = Flatten[Table[Table[{m, l}, {l, Dimensions[mflist[[m]]][[1]]}], {m, \[Mu]}], 1]; 
                       DDist = Table[Table[0, {g, Dimensions[mflist[[m]]][[1]]}], {m, \[Mu]}]; 
                       grcs = Table[Table[{}, {g, Dimensions[mflist[[m]]][[1]]}], {m, \[Mu]}]; 
                       For[g = 1, g <= \[Gamma], g++, staken = {}; 
                           For[connection = 1, connection <= d, connection++, 
                               s = RandomSample[Switch[mossyrepeat, 1, 
                               	                       Complement[glomlist, staken], 
                               	                       0, 
                                                       Select[glomlist, MemberQ[Complement[glomlist[[All, 1]], staken[[All, 1]]], #[[1]]] &], 
                                                       2, 
                                                       gl = Flatten[Table[RandomSample[Table[{m, l}, {l, Dimensions[mflist[[m]]][[1]]}], 1], {m, \[Mu]}], 1]; 
                                                       Select[gl, MemberQ[Complement[gl[[All, 1]], staken[[All, 1]]], #[[1]]] &]
                                                      ]
                                               , 1][[1]]; 
                                               AppendTo[staken, s]; 
                                               DDist[[s[[1]]]][[s[[2]]]] = DDist[[s[[1]]]][[s[[2]]]] + 1]]; 
                                               Switch[debug, 1, Print[DDist], 0, Null]; 
                                               For[connection = 1, connection <= Max[DDist], connection++, 
                                                   For[m = 1, m <= \[Mu], m++, 
                                                   	   gloms = Dimensions[DDist[[m]]][[1]]; 
                                                       For[gl = 1, gl <= gloms, gl++, 
                                                           thisgrc = If[connection <= DDist[[m]][[gl]], 
                                                           	            closestallowedGRC[grcs[[m]][[gl]], grclist, mflist[[m]][[gl]], connlist, d, dlen, omit[[m]]],
                                                           	            {-1}
                                                           	           ]; 
                                                           If[thisgrc[[1]] > 0, AppendTo[grcs[[m]][[gl]], thisgrc[[1]]], Null]; 
                                                           If[thisgrc[[1]] > 0, AppendTo[connlist[[thisgrc[[1]]]], {m, gl}], Null]; 
                                                           Switch[mossyrepeat, 1, 
                                                           	      Null, 
                                                           	      0, 
                                                                  If[thisgrc[[1]] > 0, AppendTo[omit[[m]], thisgrc[[1]]], Null]
                                                                 ]; 
                                                           If[thisgrc[[1]] == 0, 
                                                              Switch[debug, 1, 
                                                                     Print[{{m, gl}, grcs[[m]][[gl]], connection, DDist[[m]][[gl]], omit[[m]]}], 
                                                                     0,  
                                                                     Null
                                                                    ]; cnt++, 
                                                              Null
                                                             ]
                                                          ]
                                                      ]
                                                  ]; 
                                               Switch[force, 0, 
                                               	      Null, 
                                               	      1, 
                                                      noncomplete = DeleteCases[Flatten[Table[Table[
                                                      	                                           If[Dimensions[grcs[[m]][[l]]][[1]] != DDist[[m]][[l]], 
                                                      	                                           	  {m, l}, 
                                                                                                      0
                                                                                                     ], 
                                                                                                    {l, Dimensions[DDist[[m]]]}], 
                                                                                              {m, \[Mu]}],
                                                                                        1], 
                                                                                0]; 
                                                      newstruct = shuffleconns[grclist, mflist, d, dlen, noncomplete, grcs, DDist, connlist, mossyrepeat]
                                                     ]; 
                                               Switch[force, 1, 
                                               	      {newstruct[[1]], DDist, grcs, cnt, noncomplete, newstruct[[2]]},
                                                      0, 
                                                     {connlist, DDist, grcs, cnt}]
                     ]
(**************************************************************************************)
grclist2DDist[grclist_] := Flatten[Table[Table[Dimensions[grclist[[m]][[gl]]][[1]], {gl, Dimensions[grclist[[m]]][[1]]}], {m, Dimensions[grclist][[1]]}]]
(**************************************************************************************)
grclist2conns[grclist_, glompos_] := Module[{out = Table[{}, {x, Max[grclist]}]}, 
                                            Table[Table[
                                            	        AppendTo[out[[#]], {m, gl}] & /@ grclist[[m]][[gl]], 
                                            	        {gl, Dimensions[glompos[[m]]][[1]]}], 
                                                  {m, Dimensions[glompos][[1]]}]; 
                                            out]  
(**************************************************************************************)                                            
uniquemfsFROMgrclist[grclist_] := Table[Table[Dimensions[Union[grclist[[x]][[y]]]][[1]], {y, Dimensions[grclist[[x]]][[1]]}], {x, Dimensions[grclist][[1]]}]
(**************************************************************************************)  
 grcrepeats[conns_] := Module[{grc, numgrc = Dimensions[conns][[1]], d, cnt = 0}, 
                              For[grc = 1, grc <= numgrc, grc++, 
                              	  d = If[# == {}, 0, #[[1]]] &@Dimensions[conns[[grc]]]; If[d > 0, cnt = cnt + If[# != d, 1, 0] &@Dimensions[Union[conns[[grc]][[All, 1]]]][[1]]]
                              	 ]; 
                              cnt
                             ]                                            
(**************************************************************************************)  
netmossycount[conns_] := Dimensions[Union[Flatten[conns[[All, All, 1]], 1]]][[1]]
(**************************************************************************************)  
denlendist[conns_, glompos_, granpos_] := Module[{ngrans = Dimensions[conns][[1]], RAWlens}, 
                                                 RAWlens = Table[Table[EuclideanDistance[
                                                                                         glompos[[conns[[g]][[d]][[1]]]][[conns[[g]][[d]][[2]]]], granpos[[g]]], 
                                                                   	   {d, Dimensions[conns[[g]]][[1]]}], 
                                                                 {g, ngrans}]; 
                                                 RAWlens
                                                ]
 (**************************************************************************************)                                                 
meandens[bigdim_, smalldim_, Nfibres_, rdist_, xmean_, ymean_, zmean_,geom1_, geom2_, reps_] := 
         Module[{md = 0, r, gp, gd}, 
         	    For[r = 1, r <= reps, r++, 
                    gp = glopos[bigdim, Nfibres, rdist, xmean, ymean, zmean, geom1]; 
                    gd = glomdens[bigdim, smalldim, gp, geom2]; 
                    md = md + gd
                   ]; 
                md = md/reps; 
                md
               ]   
(**************************************************************************************)                
targetmf[targetglom_, biglen_, smalllen_, rdist_, glox_, gloy_, gloz_,geom1_, geom2_, reps_] := 
         Module[{numglom = 0, MAX = 50, cnt = 0, Nfibres = 5000, tol = 0.05, \[Delta] = 0.001}, 
                While[Abs[numglom - targetglom] > tol*targetglom && cnt < MAX, 
                      numglom = meandens[biglen, smalllen, Nfibres, rdist, glox, gloy, gloz, geom1, geom2, reps]; 
                      Nfibres = Nfibres + \[Delta]*(targetglom - numglom); 
                      cnt++
                     ]; 
                If[cnt >= MAX, Print["Maximum iterations reached"]]; 
                {Nfibres, numglom}
               ]
(**************************************************************************************)  
corrvsdist[shellinc_, shellmax_, nets_, rdist_, slen_, centrallen_, Nfib_, Ngrc_, glox_, gloy_, gloz_, denlen_, d_, pin_, \[Phi]_, tc_, geom1_, geom2_, mode_, d\[Mu]_: 15, d\[Sigma]_: 5] := 
           Module[{algo = 1, glompos, granpos, conns, i, A, cmat, dmat, out = Table[0, {x, nets}]}, 
           For[i = 1, i <= nets, i++, 
               glompos = glopos[slen, Nfib, rdist, glox, gloy, gloz, geom1]; 
               granpos = grcpos[slen, centrallen, Ngrc, geom2]; 
               conns = Switch[algo, 0, 
                              connections[granpos, glompos, denlen, d, mode, d\[Mu], d\[Sigma]], 
                              1, ALGconnections[denlen, glompos, granpos, d][[1]]
                             ]; 
               A = Amatrix[conns]; 
               cmat = cellcor[tc, \[Phi], pin, A]; 
               dmat = celldist[granpos]; 
               out[[i]] = Table[N[If[Dimensions[#][[1]] == 0, 0, 
               	                Mean[#]] &@ Flatten[cmat][[Flatten[Position[IntervalMemberQ[Interval[{sh - shellinc, sh}], Flatten[dmat]], True]]]
               	                ]], 
               	               {sh, shellinc, shellmax, shellinc}]
              ];
              (*Mean[out]*){out,cmat}
          ]
(**************************************************************************************)          
meangloms[nets_, rdist_, slen_, clen_, Nfib_, Ngrc_, glox_, gloy_, gloz_, denlen_, d_, geom1_, geom2_, mode_, d\[Mu]_: 15, d\[Sigma]_: 5, repeat_, force_] := 
          Module[{glompos, granpos, conns, i, gloms = 0, outside = 0, algo = 1}, 
                 For[i = 1, i <= nets, i++, 
                     glompos = Switch[outside, 0, 
                                      DeleteCases[confineglom[glopos[slen, Nfib, rdist, glox, gloy, gloz, geom1], geom2, clen, denlen], {}], 
                                      1, glopos[slen, Nfib, rdist, glox, gloy, gloz, geom1]
                                     ]; 
                     granpos = grcpos[slen, clen, Ngrc, geom2]; 
                     conns = Switch[algo, 1, 
                                    ALGconnections[denlen, glompos, granpos, d, repeat, force][[1]], 
                                    0, 
                                    connections[granpos, glompos, denlen, d, mode]
                                   ]; 
                     gloms = gloms + Dimensions[Union[Flatten[conns[[All, All, 1]]]]][[1]]
                    ]; 
                 N[gloms/nets]
                ]
(**************************************************************************************)             
Amatrix[conns_] := Module[{connectedmfs, nmfs, ngrc, A, g, Aout}, 
                           connectedmfs = Union[Flatten[conns[[All, All, 1]]]]; 
                           ngrc = Dimensions[conns][[1]]; nmfs = Dimensions[connectedmfs][[1]];
                           A = Table[Table[0, {mf, nmfs}], {grc, ngrc}]; 
                           For[g = 1, g <= ngrc, g++,
                               Table[A[[g]][[Position[connectedmfs, conns[[g]][[x]][[1]]][[1]][[1]]]] = 1, {x, Dimensions[conns[[g]]][[1]]}]
                              ]; 
                           Aout = Table[Table[If[x < nmfs && y > nmfs , A[[y - nmfs]][[x]], If[x > nmfs && y < nmfs, A[[x - nmfs]][[y]], 0]], 
                           	                  {x, nmfs + ngrc}], 
                           	            {y, nmfs + ngrc}]; 
                           {A, Aout}
                          ]
(**************************************************************************************)                         
neighbourhood[A_] := Module[{Neighbours, connlist, Nout, Nconns, n, Ngrc = Dimensions[A[[1]]][[1]]}, 
	                        connlist = Position[A[[1]], 1]; 
                            Nconns = Dimensions[connlist][[1]]; 
                            Neighbours = Table[Table[0, {x, Ngrc}], {y, Ngrc}]; 
                            For[n = 1, n <= Nconns, n++, 
                                Neighbours[[connlist[[n]][[1]]]][[connlist[[Position[connlist[[All, 2]], connlist[[n]][[2]]][[1]][[1]], 1]]]] = 1
                               ]; 
                            (*Neighbours;*) 
                            Nout = Table[Table[If[x < Ngrc && y > Ngrc , Neighbours[[y - Ngrc]][[x]], If[x > Ngrc && y < Ngrc, Neighbours[[x - Ngrc]][[y]], 0]], 
                            	               {x, 2*Ngrc}], 
                            	         {y, 2*Ngrc}]; 
                            {Neighbours, Thread[If[Thread[# >= 1], 1, 0]] & /@ (Neighbours + Transpose[Neighbours])}
                           ]
(**************************************************************************************)                           
degreedist[nets_, rdist_, slen_, clen_, Nfib_, Ngrc_, glox_, gloy_, gloz_, denlen_, d_,geom1_, geom2_, mossyrepeat_, force_] := 
           Module[{glompos, granpos, conns, i, A, dist1 = {}, dist2 = {}, outside = 0}, 
                  For[i = 1, i <= nets, i++, 
                      glompos = Switch[outside, 0, 
                                       DeleteCases[confineglom[glopos[slen, Nfib, rdist, glox, gloy, gloz, geom1], geom2, clen], {}], 
                                       1, 
                                       glopos[slen, Nfib, rdist, glox, gloy, gloz, geom1]
                                      ]; 
                      granpos = grcpos[slen, clen, Ngrc, geom2]; 
                      conns = ALGconnections[denlen, glompos, granpos, d, mossyrepeat, force][[1]]; 
                      A = Amatrix[conns]; 
                      dist1 = AppendTo[dist1, Total[Binmat[Dimensions[A[[1]]][[1]], Dimensions[A[[1]]][[2]], d]]]; 
                      dist2 = AppendTo[dist2, Total[A[[1]]]]]; 
                  {Flatten[dist1], Flatten[dist2]}]
(**************************************************************************************)                  
(*internal switch 'outside' determines whether granule cells connections are permitted to extend beyond the granule cell volume*)
anatomicalent[nets_, slen_, Nfib_, rdist_, glox_, gloy_, gloz_, geom1_, centrallen_, Ngrc_, geom2_, denlen_, d_, mode_, d\[Mu]_, d\[Sigma]_, tc_, \[Phi]_, pin_, reps_, algo_, repeats_, force_] := 
              Module[{glompos, granpos, conns, i, A,out = Table[0, {x, nets}], outside = 0}, 
                     For[i = 1, i <= nets, i++, 
                         Switch[outside, 0, 
                         	    glompos = DeleteCases[confineglom[glopos[slen, Nfib, rdist, glox, gloy, gloz, geom1], geom2, centrallen, denlen], {}], 
                         	    1, glompos = glopos[slen, Nfib, rdist, glox, gloy, gloz, geom1]
                         	   ]; 
                         granpos = grcpos[slen, centrallen, Ngrc, geom2]; 
                         conns = Switch[algo, 2, 
                         	            connections2[granpos, glompos, denlen, d], 
                                        1, 
                                        ALGconnections[denlen, glompos, granpos, d, repeats, force], 
                                        0, 
                                        connections[granpos, glompos, denlen, d, mode, d\[Mu], d\[Sigma]]
                                       ]; 
                         A = Amatrix[conns[[1]]][[1]]; 
                         out[[i]] = FindentNUM[tc, Dimensions[A][[2]], d, \[Phi], Dimensions[A][[1]], pin, reps, 1, Transpose[A]]
                        ]; 
                        Mean[out]
                       ]
(**************************************************************************************)
tissueprof[rdist_, bigdim_, cendim_, nmfs_, dlen_, d_, glox_, gloy_, gloz_, geom1_, geom2_, mode_, force_] := 
           Module[{glomraw, glompos, granpos, conns},
                   glomraw = glopos[bigdim, nmfs, rdist, glox, gloy, gloz, geom1];
                   glompos = DeleteCases[confineglom[glomraw, geom2, cendim], {}];
                   granpos = grcpos[bigdim, cendim, grcinvol[cendim, geom2], geom2];
                   conns = ALGconnections[dlen, glompos, granpos, d, mode, force]; 
                   {Dimensions[glompos][[1]], Dimensions[Union[Flatten[conns[[1]][[All, All, 1]]]]][[1]], 
                    grcrepeats[conns[[1]]], uniqueconns[conns], Dimensions[Union[Flatten[conns[[1]], 1]]][[1]], 
                    Dimensions[Cases[Table[Dimensions[conns[[1]][[x]]][[1]], {x, Dimensions[conns[[1]]][[1]]}] - d, 0]][[1]], 
                    Dimensions[Cases[Table[Dimensions[Union[conns[[1]][[x]]]][[1]], {x, Dimensions[conns[[1]]][[1]]}] - d, 0]][[1]] - 
                    Dimensions[Cases[Table[Dimensions[conns[[1]][[x]]][[1]], {x, Dimensions[conns[[1]]][[1]]}] - d, 0]][[1]], 
                    BinCounts[Table[Dimensions[Union[conns[[1]][[x]][[All, 1]]]][[1]], {x, Dimensions[conns[[1]]][[1]]}]]}
                 ]
(**************************************************************************************)
meanconpro[rdist_, bigdim_, cendim_, nmfs_, dlen_, d_, glox_, gloy_, gloz_, geom1_, geom2_, reps_, mode_, force_] := 
           Module[{glomraw, glompos, granpos, conns, gr = 0, uc = 0, rep},
                  For[rep = 1, rep <= reps, rep++, 
                  	  glomraw = glopos[bigdim, nmfs, rdist, glox, gloy, gloz, geom1];
                      glompos = DeleteCases[confineglom[glomraw, geom2, cendim], {}];
                      granpos = grcpos[bigdim, cendim, grcinvol[cendim, geom2], geom2];
                      conns = ALGconnections[dlen, glompos, granpos, d, mode, force]; 
                      gr = gr + grcrepeats[conns[[1]]]; 
                      uc = uc + Mean[Table[Dimensions[Union[conns[[1]][[x]][[All, 1]]]][[1]], {x, Dimensions[conns[[1]]][[1]]}]]
                     ]; 
                  {gr/reps, uc/reps}
                 ]
(**************************************************************************************)          
End[]
EndPackage[]

