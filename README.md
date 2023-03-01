# querysa 
This program is for the second half part of the course CMSC701: Computational Genomics Programming Assignment1.

Implement a program that reads in the serialized data structures and then performs a series of queries against the suffix array.

This program will take as input 4 arguments, your serialized suffix array, a `FASTA` file containing queries, a `query mode` parameter, and an `output` file name. It will perform query in the SA and report the results in the `output` file specified in the format specified below.

## Input

- `index` - the path to the binary file containing your serialized suffix array (as written by `buildsa`).

- `queries` - the path to an input file in `FASTA` format containing a set of records. 

- `query mode` - this argument should be one of two strings; either `naive` or `simpaccel`. If the string is `naive`, the program will perform your queries using the `naive` binary search algorithm. If the string is `simpaccel`, the program can perform your queries using the “simple accelerant” algorithm we covered in class.

- `output` - the name to use for the resulting output.

## Output


- `output` - the output file of your program. This file should contain the results of your queries in the following format. Each line should contain a **tab-separated** list containing the following information:
```
query_name, k, hit_1, hit_2, hit_k
```

## Example

To run this program on Linux, download the whole `querysa` folder, go to the `querysa` directory and make sure the parameter given by this order: `index`, `queries`, `query_mode`, `output`. 

Here is an example:

```{bash}
$ cargo run output_ecoli_12mer.bin query.fa naive result_ecoli_12mer.txt
```