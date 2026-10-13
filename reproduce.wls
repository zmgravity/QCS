root=DirectoryName[$InputFileName];
Get[FileNameJoin[{root,"Kernel","QCS.wl"}]];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->1];
args=If[Length[$ScriptCommandLine]>0,Rest[$ScriptCommandLine],
 With[{position=FirstPosition[$CommandLine,"-script"]},If[MissingQ[position],{},Drop[$CommandLine,First[position]+1]]]];
workflow=If[Length[args]>0,First[args],"check"];
result=QCSRun[workflow];
If[result===$Failed,Quit[1]];
Print["Finished workflow: ",workflow,". Results: ",FileNameJoin[{QCSRoot,"results"}]];
Quit[0];
