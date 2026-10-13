root=DirectoryName[$InputFileName];
lines=Select[Import[FileNameJoin[{root,"MANIFEST.sha256"}],"Lines"],StringLength[#]>0&];
failed=Select[lines,Function[line,With[{hash=StringTake[line,64],relative=StringDrop[line,66]},
 !FileExistsQ[FileNameJoin[Prepend[StringSplit[relative,"/"],root]]]||
 IntegerString[FileHash[FileNameJoin[Prepend[StringSplit[relative,"/"],root]],"SHA256"],16,64]=!=hash]]];
Print[Length[lines]-Length[failed]," / ",Length[lines]," files match the manifest."];
If[Length[failed]>0,Print[failed];Quit[1],Quit[0]];
