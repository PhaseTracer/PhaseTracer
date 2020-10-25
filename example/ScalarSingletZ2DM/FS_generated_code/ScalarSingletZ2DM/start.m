AppendTo[$Path, FileNameJoin[{Directory[], "meta"}]];

Needs["SARAH`"];
Needs["FlexibleSUSY`", FileNameJoin[{Directory[], "meta", "FlexibleSUSY.m"}]];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = FileNameJoin[{workingDirectory, "Output"}];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{workingDirectory, "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};

Print["Current working directory: ", workingDirectory];
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];
Print["Starting model: ", ToString[HoldForm[Start["ScalarSingletZ2DM"]]]];

Start["ScalarSingletZ2DM"];

MakeFlexibleSUSY[InputFile -> FileNameJoin[{"models/ScalarSingletZ2DM","FlexibleSUSY.m"}],
                 OutputDirectory -> "models/ScalarSingletZ2DM",
                 DebugOutput -> False];
