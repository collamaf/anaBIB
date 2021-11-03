# Software package to perform BIB analysis
*Last update: 28-10-2021 by collamaf*

----------

The package is composed by 2 python programs:

`ReadBIBfiles`:
- Reads the `DUMP` files produced by FLUKA and produces a single text file with all the relevant info needed to perform BIB analysis


`AnaBIB`
- Performs the actual analysis

## How to run

- To read fluka `DUMP` files:
```
python --dataPath "/path/to/DUMPfiles/with*" --outPath "/path/to/outputfile" ReadBIBFiles.py
```
[please note the "" if you want to use wildcards (*...)]


- To analyse BIB:
```
python anaBib.py --runName "MyRun" --fileList File1 File2 --labelList Label1 Label2 --ele
```

The program reads the files (1 or more), and creates several `.png ` files of relevant plots (comparing variables per particle and/or per run). 

The `--ele` flag (default False) requests plots based on the `XX_ele` output file with info regarding parent electron. 

Png files of the plots are produced in the same folder.

