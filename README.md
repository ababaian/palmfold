```
  _____      _            __      _     _ 
 |  __ \    | |          / _|    | |   | |
 | |__) |_ _| |_ __ ___ | |_ ___ | | __| |
 |  ___/ _` | | '_ ` _ \|  _/ _ \| |/ _` |
 | |  | (_| | | | | | | | || (_) | | (_| |
 |_|   \__,_|_|_| |_| |_|_| \___/|_|\__,_|
 =========================================
 v221108
```

Structural classifier of RNA viral RNA-dependent RNA polymerase (RdRp). Given a protein structure file (.pdb), `palmfold` performs structural alignments against a constellation of RdRp and non-RdRp polymerase (XdXp) structures and scores the matches. It will evaluate if the input structure is predicted to be RdRp, independent of sequence information.

## Usage

```
python3 palmfold.py -i ./pdb -p ./pol -o ./test

arguments:

  -h, --help`
  						show this help message and exit
  						
  --inputpath INPUTPATH, -i INPUTPATH
                        Directory containing "*.pdb[.gz]" files to be classified. [./pdb]

  --palmprints PALMPRINTS, -p PALMPRINTS
                        path to "./palmfold/pol/" containing palmprint pdb models. [./pol]

  --outpath OUTPATH, -o OUTPATH
                        Directory for writing output files. [./out]

  --threshold THRESHOLD, -t THRESHOLD
                        Polymerase+ classification threshold for TMAlign. [0.5]

  --verbose, -v
    					Print INFO messages, else prints WARN and ERROR only
```


## Install

### Dependencies

- `Python3` & `Bio` library : `pip3 install "Bio"`
- [`TMalign`](https://zhanggroup.org//TM-align/) ([binary link](https://zhanggroup.org//TM-align/TMalign.gz)), available in `$PATH`

```
# Clone respository
git clone https://github.com/ababaian/palmfold.git

cd palmfold

# Unit test for palmfold
python3 palmfold.py -i ./pdb -p ./pol -o ./test
```

