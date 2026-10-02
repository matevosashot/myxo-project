(* ::Package:: *)

(* RemoteExport: Export that works with a remote kernel.

   With a remote kernel, Export hands PDF/EPS/SVG/... rendering to the
   front end, which then tries to write the file on the front end's
   machine at the kernel's path, and fails silently. RemoteExport runs
   Export in a subkernel on the kernel's machine instead; that subkernel
   has no front end attached, so it renders with its own headless front
   end and writes the file locally.

   Usage:
     Needs["RemoteExport`"]
     RemoteExport["a.pdf", Plot[Sin[x], {x, 0, 10}]]
*)

BeginPackage["RemoteExport`"]

RemoteExport::usage =
  "RemoteExport[file, expr] exports expr to file, rendering it in a subkernel on the kernel's machine.\n" <>
  "RemoteExport[file, expr, format, opts] passes format and options to Export.";

RemoteExport::nokernel = "Could not launch the export subkernel.";

Begin["`Private`"]

$exportKernel = None;

(* (re)launch the subkernel if it is not running, e.g. after CloseKernels[] *)
exportKernel[] := (
  If[!MemberQ[Kernels[], $exportKernel],
    $exportKernel = Replace[LaunchKernels[1], {{k_} :> k, _ :> None}]];
  $exportKernel)

RemoteExport[file_String, expr_, args___] :=
  Module[{k = exportKernel[], out = ExpandFileName[file], res},
    If[k === None,
      Message[RemoteExport::nokernel]; $Failed,
      (* absolute path, since the subkernel's Directory[] may differ *)
      res = With[{o = out, e = expr, a = {args}},
        ParallelEvaluate[Export[o, e, Sequence @@ a], k]];
      If[StringQ[res] && FileExistsQ[out], out, $Failed]]]

(* launch the subkernel when the package is loaded *)
exportKernel[];

End[]

EndPackage[]
