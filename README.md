# vspc_graph

Find the set of all rings (cycles) in an undirected graph. The algorithm was inspired by the work presented here, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2765087/, which will be properly cited soon. However, we've had to take our own liberties with adapting their algorithm since it currently seems like we aren't necessarily looking for an SSSR (smallest set of smallest rings). The code still uses the term SSSR because it was how it all began. This will be modified soon.

## Requirements

You will need the following in order to build and run this software:
* CMake 3.10.2 or higher
* C++17-compliant compiler
* Python3 interpreter 

The deatils on reaching these requirements is omitted here since a few Google searches should do the trick. 

## How to build

First, either clone this repository (or your own fork of it), or download and unzip the latest release from the release page. Then, run the following commands (replacing anything surrounded by the angled brackets, including the brackets, with the appropriate paths):

```
cd <path/to/this/local/repo>
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=<path/to/C++/compiler>
make install
```

This will build an executable from the C++ code, called `main`, and install it in the `bin` directory of this repository. The `bin` directory is also where you'll find `sssr.py`, the Python-wrapper script for the C++ executable that is designed to help users interface with the executable. **Important**: while `main` and `sssr.py` need not live in the `bin` directory of this respository in order for everything to work, they must live in the same directory. 

## How to use

Because of the mediocore installation process that is currently implemented, we highly suggest adding a command-line alias to make it much more convenient to call the Python script from anywhere, e.g.
```
alias sssr='<path/to/this/local/repo/bin/sssr.py'
```
A Google search shouold be sufficient for helping you get this set up. 

To use this software, run something like the following:
```
sssr path/to/data/graph_file.csv
``` 

This will process the graph data in `path/to/data/graph_file.csv` and produce a CSV file containing cycle data called `path/to/data/cycles_file.csv`. A few things to note:
* Files that do not have the csv extension are automatically ignored.
* Files that do not start with "graph" are automatically ignored.
* Files that begin with "cycles" are automatically ignored. The reason for this one is because the software creates the names of the output files by replacing "graph" of the input file name with "cycles. Hence, it sees files that begin with "cycles" as one of its own output files.

The Python script supports the input of multiple files at once, even via the use of wildcards. Furthermore, multiple files, by default, are processed in parallel! For example,
```
sssr path/to/data/*.csv path/to/more/data/*.csv
```
will attempt and process all of the CSV files in `path/to/data` and `path/to/more/data` according to the three criteria above. Note that a slightly shorter version of the above command would also work, thanks to the first criterion above:
```
sssr path/to/data/* path/to/more/data/*
```
Depending on how many files are in the two directories combined, a possible downside to this version is that the python script might have to do a bit more work to filter out the non-CSV files. 

As mentioned above, given a valid input file name, the software will automatically generate the output file name by replace "graph" with "cycles," e.g. an input file name of `graph_test.csv` will produce an output file called `cycles_test.csv` in the same directory. Thus, the software is able to check whether or not any of the input files have already been processed (by checking for the existence of its associated output file). If so, it will ask whether or not you would like to overwrite the output file. Enter either `y` or `n` for each of the cases. This feature allows you to not have to worry about being overly specific when specifying input files, e.g. purposely ommitting graph files that you don't need to reprocess. 

In future releases, we plan to add more features to the Python script. The following are just a few examples:
* The option to specify how many threads you want to dedicate to this entire process.
* The option to automatically overwrite all existing output files with newly generated results.

## Future work

We will work on a better process of installation once it's confirmed that a Python version of Lindsay's code is in the works. Currently, this piece of software is isolated from the rest of the code, so there doesn't seem to be an immediate need to develop a proper installation process.

It was mentioned earlier that while we did follow an existing algorithm to develop the C++ portion of the code, we had to make some modifications to fit the needs of Lindsay's project. This calls for some properly written documentation for what our version of the algorithm is doing. 

