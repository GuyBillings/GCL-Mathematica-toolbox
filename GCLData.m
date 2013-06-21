(*Mathematica Package to load Billings and Silver 2013 data relating the cerebellar granule cell layer*)
(* Created by the Wolfram Workbench Jun 19, 2013 *)

BeginPackage["GCLData`"]
Needs["BinaryPartitionNetworks`"]
(* Exported symbols added here with SymbolName::usage *) 
DataInfo::usage="Info[]: Print information about the data in this package"
LoadGCLdataset::usage="LoadGCLdataset[path_,dataset_]: Load precomputed GCL data from the root path 'path' (see readme accompanying data for instructions).\
specify dataset to load with 'dataset'. Default to load is GCLpaperdataset1."
(**************************************************************************************)
Begin["`Private`"]
(* Implementation of the package *)
GCLpaperdataset1="Theory data set for 1 billion patterns: Results of entropy calculation method for a network representing an 80 micron diameter ball of GCL.\ 
Network contains 176 MF inputs and 509 GC outputs. Network encodes 1 billion samples (i.e. 1 billion trials of the stochastic MF variables are encoded).\
Dataset is for this encoding over MF activity levels from 0.005 to 0.995 at 0.005 increments and 1 to 20 connections per GC for all threshold levels."
GCLpaperdataset2="Theory data set for 4000 patterns: Results of entropy calculation method for a network representing an 80 micron diameter ball of GCL.\ 
Network contains 176 MF inputs and 509 GC outputs. Network encodes 4000 samples (i.e. 4000 trials of the stochastic MF variables are encoded).\
Dataset is for this encoding over MF activity levels from 0.005 to 0.995 at 0.005 increments and 1 to 20 connections per GC for all threshold levels."
GCLpaperdataset3="Anatomical model and graph data set for 4000 patterns: Results of simulation of the anatomical/graph model network representing an 80 micron diameter ball of GCL.\ 
Network contains 176 MF inputs and 509 GC outputs on average. Network encodes 4000 samples (i.e. 4000 of the stochastic MF variables are encoded).\
Dataset is for this encoding over MF activity levels from 0.005 to 0.995 at 0.005 increments and 1 to 20 connections per GC for all threshold levels. An \
average was taken over 10 reptetitions for each parameter configuration.\
Loads array having format: For each p(MF) and d value there is {matrix_model_constrained[[\[Phi]]], anatomical_model _constrained[[\[Phi]]], \
matrix_model _unconstrained[[\[Phi]]], anatomical_model _unconstrained[[\[Phi]]]}, where constrained is for models having connectivity consraints and\
unconstrained is for models with out the constrainsts (see sumplimentary material for Billings Silver 2013)."
GCLpaperdataset4="Bipartite graph model data set for 4000 patterns with MF correlations as described in Billings & Silver 2013: \
Results of simulation of the Bipartite graph model network representing an 80 micron diameter ball of GCL.\ 
Network contains 176 MF inputs and 509 GC outputs on average. Network encodes 4000 samples (i.e. 4000 of the stochastic MF variables are encoded).\
Dataset is for this encoding over MF activity levels from 0.005 to 0.995 at 0.005 increments and 1 to 20 connections per GC for all threshold levels. An \
average was taken over 10 reptetitions for each parameter configuration."

loadmodeldata[path_, base_, suf_, R_, repeat_, d_] := 
             Module[{tmp, 
             	     p, 
             	     debug = 0, 
             	     fullpath, 
             	     out, 
             	     fns}, 
                     fns = FileNames[path <> ToString[R] <> "_" <> ToString[repeat] <> "_" <> ToString[d] <> "/" <> "*" <> suf]; 
  					 out = Table[0, {x, Dimensions[fns][[1]]}]; 
  					 For[p = 1, p <= Dimensions[fns][[1]], p++, fullpath = fns[[p]]; 
  					 	             Switch[debug, 0, Null, 1, Print[fullpath]]; 
   									 tmp = Partition[Flatten[ToExpression[StringSplit[ReadList[fullpath, "String"], "\t"]]], d]; out[[p]] = tmp]; out]
   									 
loadmodeldata2[path_, base_, suf_, R_, repeat_, d_] := 
 			 Module[{tmp, 
 			 	     p, 
 			 	     debug = 0, 
 			 	     fullpath, 
 			 	     out, 
 			 	     fns}, 
  					 fns = FileNames[path <> "80_" <> ToString[d] <> "_1/" <> base <> ToString[R] <> "_" <> "*" <> "_" <> ToString[d] <> "_" <> ToString[repeat] <> suf];
  					 out = Table[0, {x, Dimensions[fns][[1]]}]; 
                     For[p = 1, p <= Dimensions[fns][[1]], p++, fullpath = fns[[p]]; 
                                                           Switch[debug, 0, Null, 1, Print[fullpath]]; 
   														   tmp = Partition[
      													   Flatten[ToExpression[StringSplit[ReadList[fullpath, "String"], "\t"]]], d]; out[[p]] = tmp]; out]   									 
(**************************************************************************************)
DataInfo[dataset_:1]:=
Module[{},
   Switch[dataset,1,
   	              Module[{},Print["GCLpaperdataset1:"];
	              Print[GCLpaperdataset1];
				  Print["(**************************************************************************************)"]],
				  2,
				  Module[{},Print["GCLpaperdataset2:"];
	              Print[GCLpaperdataset2];
				  Print["(**************************************************************************************)"]],
				  3,
				  Module[{},Print["GCLpaperdataset3:"];
	              Print[GCLpaperdataset3];
				  Print["(**************************************************************************************)"]],
				  4,
				  Module[{},Print["GCLpaperdataset4:"];
	              Print[GCLpaperdataset4];
				  Print["(**************************************************************************************)"]]
         ]]
(**************************************************************************************)
LoadGCLdataset[path_,dataset_:"GCLpaperdataset1"]:=Switch[dataset,
	  														"GCLpaperdataset1",
	  													    Loaddump[path<>"176_509_entropy.dat"],
	  													    "GCLpaperdataset2",
	  													    Loaddump[path<>"176_509_entropy_4000events.dat"],
	  													    "GCLpaperdataset3",
	  													    Table[Table[loadmodeldata[path,"BIGNET_2012-06-29.11-51-38_",".dat",80,repeat,d], {d, 20}], {repeat, 10}],
	  													    "GCLpaperdataset4",
	  													    Table[Table[loadmodeldata2[path,"correlations_2013-05-24.17-57-46_",".dat",80,repeat,d], {d, 20}], {repeat, 1,1,1}]
														 ]
(**************************************************************************************)
End[]
EndPackage[]

