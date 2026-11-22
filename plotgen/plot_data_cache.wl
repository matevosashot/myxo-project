(* ::Package:: *)

(* plot_data_cache.wl -- plot data kept next to the figure scripts (2026-11-22).

   The simulation output in <repoRoot>/data is untracked and large. A figure script reads from
   it only what it plots, stores that in a cache file in plotgen, and runs from the cache when
   the source is not there:

     data = plotDataCache[cacheFile, {source1, source2, ...}, extract]

   extract[] reads the sources and returns the plot data (any expression). Then
   - all sources exist: each is fingerprinted by its name, size and SHA-256 hash. If the cache
     holds the same fingerprints, its data are returned without reading the sources ("source
     unchanged"). Otherwise extract[] runs and the cache is rewritten ("source new", "source
     changed", or "other source" when the cache was made from different files).
   - a source is missing: the cache is used if it was made from sources of the same names
     ("source not found ... cache of <date>"); if not, the script stops with exit code 1 and
     names the missing file, since there is nothing to plot.
   Source names are relative to the repository root (the parent of plotgen), so a cache made on
   one machine is valid on another. The cache is a WXF file (Export[..., "WXF"]: exact numbers,
   readable by any Mathematica version, no front end needed) holding
     <|"sources" -> {<|"file", "bytes", "sha256"|>, ...}, "created" -> date, "data" -> data|>.
   To look into one: Import["plotdata_<tag>.wxf"].
   Usage:  Get[FileNameJoin[{<repoRoot>, "plotgen", "plot_data_cache.wl"}]] *)

$plotDataCacheRoot = ParentDirectory[DirectoryName[$InputFileName]];

plotDataSourceName[f_String] := Module[{full = ExpandFileName[f],
    root = ExpandFileName[$plotDataCacheRoot] <> $PathnameSeparator},
   If[StringStartsQ[full, root],
    StringRiffle[FileNameSplit[StringDrop[full, StringLength[root]]], "/"], full]];

plotDataFingerprint[f_String] := <|"file" -> plotDataSourceName[f], "bytes" -> FileByteCount[f],
   "sha256" -> FileHash[f, "SHA256", All, "HexString"]|>;

plotDataCache[cacheFile_String, sources : {__String}, extract_] := Module[
  {cache = None, names, missing, fp, status, data, cacheName = FileNameTake[cacheFile]},
  names = plotDataSourceName /@ sources;
  If[FileExistsQ[cacheFile], cache = Quiet@Check[Import[cacheFile, "WXF"], None]];
  If[! AssociationQ[cache] || ! SubsetQ[Keys[cache], {"sources", "created", "data"}], cache = None];
  missing = Select[sources, ! FileExistsQ[#] &];
  Which[
   (* the sources are here: use the cache only if they are the ones it was made from *)
   missing === {},
   fp = plotDataFingerprint /@ sources;
   If[cache =!= None && cache["sources"] === fp,
    Print["plot data: ", cacheName, " (source unchanged: ", StringRiffle[names, ", "], ")"];
    cache["data"],
    status = Which[cache === None, "source new",
      Lookup[cache["sources"], "file"] === names, "source changed",
      True, "other source (the cache was made from " <>
       StringRiffle[Lookup[cache["sources"], "file"], ", "] <> ")"];
    Print["plot data: ", status, ": reading ", StringRiffle[names, ", "], " and writing ", cacheFile];
    data = extract[];
    Export[cacheFile, <|"sources" -> fp, "created" -> DateString["ISODateTime"], "data" -> data|>,
     "WXF"];
    data],
   (* no source, but a cache of it *)
   cache =!= None && Lookup[cache["sources"], "file"] === names,
   Print["plot data: ", cacheName, " (source not found: ",
    StringRiffle[plotDataSourceName /@ missing, ", "], "; cache of ", cache["created"], ")"];
   cache["data"],
   (* nothing to plot *)
   True,
   Print["*** plot data: source not found: ", StringRiffle[missing, ", "], "\n    and ",
    If[cache === None, "no cache " <> cacheFile,
     "the cache " <> cacheFile <> " was made from " <> StringRiffle[Lookup[cache["sources"], "file"], ", "]],
    ".\n    Provide the simulation file."];
   Exit[1]]];
