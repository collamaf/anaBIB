# Software package to perform BIB analysis
*Last update: 26-07-2021 by collamaf*

----------

The package is composed by 2 python programs:

`ReadBIBfiles`:
- Reads the `DUMP` files produced by FLUKA and produces a single text file with all the relevant info needed to perform BIB analysis


`AnaBIB`
- Performs the actual analysis

## How to run

- To read fluka `DUMP` files:
```
python --dataPath /path/to/DUMPfiles/with* --outPath /path/to/outputfile ReadBIBFiles.py
```

- To analyse BIB:
```
python anaBib.py
```


This program looks for both 1.5 and 3TeV Monte Carlo results in `DigFiles/Part1.5TeV.dump` and `DigFiles/Part3TeV.dump` respectively.
Png files of the plots are produced in the same folder.
