<img height=40 src='http://my5C.umassmed.edu/images/3DG.png' title='3D-Genome' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/dekkerlabbioinformatics.gif' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/umasslogo.gif' />

# tab2hdf

convert tsv (tab seperated value) text (matrix) files into hdf5 (h5) files

```
tab2hdf/
    tab2hdf.py - convert tsv matrix file into hdf5 file
```

## Installation
tab2hdf requires numpy/h5py

You can install the dependencies with:
```
for req in $(cat requirements.txt); do pip install $req; done
```

Download the project.
```
wget -O tab2hdf.zip https://github.com/blajoie/tab2hdf/archive/master.zip
```
Or clone the git project
```
[ssh] - git clone git@github.com:blajoie/tab2hdf.git
[https] - git clone https://github.com/blajoie/tab2hdf.git
```

Unzip the master:
```
unzip tab2hdf.zip
cd tab2hdf/
```

After installing dependencies, you should be free to use tab2hdf.py:
e.g.
```
$ python scripts/tab2hdf.py
```

## Full Documentation

See the [tab2hdf Wiki](https://github.com/blajoie/tab2hdf/wiki) for full documentation of Hi-C HDF5 file format.
<br>

See the [HDF5 spec Wiki](https://github.com/blajoie/tab2hdf/wiki/H5-Spec) for focumentation of the dekkerlab Hi-C HDF5 file format. (h5 dataset/attributes)

Download/Clone the [tab2hdf.py](https://github.com/blajoie/tab2hdf) HDF5 helper script (python).
<br>
numpy/scipy/h5py required. Python 2.7+

See the [Usage Wiki](https://github.com/blajoie/tab2hdf#usage</a>) for help running the tab2hdf.py scripy.

See the [Wiki](https://github.com/blajoie/tab2hdf/wiki) for full documentation, examples, operational details and other information.

## Communication

- [Bryan Lajoie](https://github.com/blajoie)
- [Noam Kaplan](https://github.com/NoamKaplan)
- Twitter: [@my5C](https://twitter.com/my5C)

## What does it do?

tab2hdf can transform a tsv matrix file into a hdf5 matrix file (for random access)

## Usage

```

$ python scripts/tab2hdf.py  --help

usage: tab2hdf.py [-h] -i INFILE [-v] [-o OUTFILE] [-b BLOCKSIZE] [--version]

Extract c-row_stripe from HDF5 file into TXT (matrix.gz)

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        interaction matrix hdf5 file (default: None)
  -v, --verbose         Increase verbosity (specify multiple times for more)
                        (default: None)
  -o OUTFILE, --output OUTFILE
                        interaction matrix output file (default: None)
  -b BLOCKSIZE, --blocksize BLOCKSIZE
                        block size of HDF5 file (default: 32)
  --version             show program's version number and exit

```

## Usage Examples

## Change Log

## Bugs and Feedback

For bugs, questions and discussions please use the [Github Issues](https://github.com/blajoie/tab2hdf/issues).

## LICENSE

Licensed under the Apache License, Version 2.0 (the 'License');
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an 'AS IS' BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

