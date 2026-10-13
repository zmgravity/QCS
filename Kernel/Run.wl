QCSRun[name_String]:=Switch[name,
 "check",QCSCheckData[],"pair",QCSReconstructPair[QCSPreparePair[]],
 "pair-from-operators",QCSReconstructPair[QCSPreparePair[],True],
 "operators",Table[With[{k=QCSKernels[n,True]},Export[path["results/kernels_n"<>ToString[n]<>".wxf"],k,"WXF"]],{n,{32,48,64}}],
 "gr",QCSGR[],"flux",QCSFlux[],"wavepacket",QCSWavepacket[],"two-pole",QCSTwoPole[],"causal",QCSCausal[],"kerr",QCSKerr[],
 "rwz",QCSRWZ[],"scattering",QCSScattering[],"pole-scattering",QCSPoleScattering[],"figures",QCSFigures[],"new-figures",QCSFigures[True],
 "all",Table[Print["Running ",item];QCSRun[item],{item,{"check","pair","gr","flux","wavepacket","two-pole","causal","kerr","rwz","scattering","pole-scattering","figures"}}],
 _,Message[QCSRun::name,name];$Failed];
QCSRun::name="Unknown workflow `1`. See README.md for workflow names.";
