---
title: "Essential 'omics (bash / linux)"
layout: single
category: page
permalink: /Misc

toc: true
toc_sticky: true

header:
  image: /assets/images/home_banner.svg
---

# Basic reference for useful 'omics-related command line

A central reference guide of some essential bash programs and their microbiology-relevant uses can be really helpful for those first getting their feet wet with sequence data. This is all about lowering barriers to working with 'omics datasets by making an easy cheat-sheet of common programs and their usage. 

Many different intro-to-coding and intro-to-bioinformatics tutorials online, some of which are incredibly detailed. This isn't meant to be that, just a reference and usage guide. Some tutorials I like are linked at the bottom.

## Potentially confusing vocab, symbols, and syntax
* special words: 'directory' is the same as 'folder'
* special characters:
  * `\` : Means interpret next character literally if it's a special character, or give a special meaning to next character if it's not special yet.
    * A file with spaces, like `my file.txt` would need to be referred to as `my\ file.txt`  
  * `\t` : tab character
  * `\n` : newline character
  * `^` : used in [regular expressions](http://web.mit.edu/hackl/www/lab/turkshop/slides/regex-cheatsheet.pdf) to signify the start of a line, with tools like `sed` and `grep`
  * `$` : End of line marker, it is the counterpart of `^`
  * `.` : 'here' in the context of a filepath, or wildcard for 'any character' in a regular expression

## Navigation
* `pwd` : **p**rint **w**orking **d**irectory (print the path to where you are)
* `cd` : **c**hange **d**irectory (move to a different directory)
  * Use it like `cd DESTINATION` where DESTINATION is the path to where you want to go 
  * There are two kinds of ways to tell how to get to DESTINATION, absolute and relative:
    * Absolute starts with a `/` character, signifying the path starts at the 'root' of the file structure. Like `cd /Users/rosalindfranklin/Documents`
    * Relative is relative to current location, so `cd my_folder` would go to `my_folder` inside the current directory (if it exists, otherwise, error!) `..` is a special character that means move 'up' (towards the root) from where you are now, like `cd ../my_folder`
* `ls` : **l**i**s**t contents of the current directory
  * `ls --color=none` - Turn off colors
  * `ls PATH` - List the contents of directory/file specified at PATH   
* `*` : asterisks in file paths stands for any character, any number of times.
  * use with `ls` and other commands to specify multiple files, like `ls *txt` to list all the files in the current directory ending with `*txt`, but exclude the others
  
## Printing file contents to the screen
* `head FILE` : print the **first** 10 lines of *FILE*, great way to peek at a big file
* `tail FILE` : print the **last** 10 lines of *FILE*
* `more FILE` : print potentially all the lines of a file, one full terminal window at a time
  * Once `more` is run, pressing the space key will advance to the next page
  *  `q` quits the printing, returning to the command prompt
* `cat FILE1 [FILE2 FILE3]` - prints the whole file, very useful for con**cat**enating multiple files
  * `cat` can take multiple files, like `cat my_file.txt my_other_file.txt` and it will print the first one, then the second one, etc.  

## Writing and modifying files / directories
* `>` : redirects the output to a file instead of printing to the screen
  * `cat my_file.txt my_other_file.txt > combined.txt` concatenates two files, top of 2nd after bottom of 1st, into a file called `combined.txt`
  * Will overwrite an existing file of the same name
* `>>` : Just like `>` except it appends to the bottom of specified file, if it already exists, or creates it if it doesn't exist already. 
* `|` : the pipe character chains commands - the output of first command doesn't get printed but is passed to next command.
  * `cat FILE | head` first prints the entire contents of FILE using `cat`, but instead of printing it, the output gets passed to `head` which prints the first 10 lines of that.
* `mkdir NAME` : make directory called NAME
* `mv SOMETHING DESTINATION` : move a file or directory at SOMETHING to DESTINATION.
  * Can rename a file or directory's like `mv file.txt newname.txt` instead of moving it into a new directory
* `cp SOMETHING DESTINATION` : make a copy of file SOMETHING at DESTINATION
* `rm FILE [FILE2 FILE3]` : removes (deletes) FILE. **BE CAREFUL -  THIS IS IRREVERSIBLE**. Ensure you have filenames typed carefully and no unexpected asterisks or spaces to not accidentally delete everything on your computer. It's happened...
  * `rm -rf DIRECTORY` : the `-rf` tells it to delete directories and their contents too

## Searching and manipulating text
### `grep`
* `grep "SUBSTRING" FILE` : the command line version of find, searches for SUBSTRING in FILE and returns only lines that contain SUBSTRING. Some useful specific applications:
* `grep -c ">" *fasta` : `-c` means 'count' instead of print, so this counts the number of sequences in each fasta file in the current directory
* `grep -A 1 ">my_sequence" many_sequences.fasta` : `-A N` returns `N` lines after the match, so if `many_sequences.fasta` is nicely formatted, this will return the sequence named `my_sequence`
  * Like `-A`, `-B N` returns N lines before the match, and `-C N` returns N lines of context (before and after the match)
* `grep -v "metazoan" organisms.txt` : `-v` means invert, so this returns all lines that DON'T contain "metazoan"
* `grep "THIS" FILE | grep "THAT"` : Chaining greps with a pipe is the equivalent of a logical AND; this returns lines containing **both** THIS and THAT
 * `grep "THIS\|THAT" FILE` : `|` is a logical OR, but in grep it needs a `/` to signify that it is not a literal "\|" character to match, so this command returns lines with THIS or THAT (or both)

### `awk`
Think Microsoft Excel of the command line, very useful for tables and columnar data. Some examples:
* `awk -F"\t" '{print $0}' FILE` : The `-F` flag lets you specify what the column separator is (in this case, tabs), then in the single quotes and curly braces is the program awk runs. `print` just means print, and `$0` means all columns
* `awk -F"\t" '{print $1}' FILE` means print only the first column of FILE, if columns are tab-separated
* `awk -F"\t" '{print $1 "_" $2}' FILE` parses FILE as tab-separated, takes only the first two columns, and pastes them together with an underscore.
* `awk -F" " '{print $1}' FASTA` is a great way to clean up a fasta if there is junk info after whitespaces in the defline, assuming the sequence doesn't have spaces (check first with `head`)
* Built-in awk variables: 
  * `NR` is rownumber 
  * `NF` is number of fields in a given row
  * `awk '{print NF}' FILE` prints the number of fields (columns) in each line (useful for identifying misbehaving rows if you get an error related to having an irregular number of columns)

### Others
* `tr ' ' '_' < FILE` : **tr**anslates space characters (first argument) to underscores (second argument). The `<` is required for `tr` to read in the file, otherwise `cat FILE | tr ' ' '_'`
* `sed "s/SEARCH/REPLACE/g" FILE` : replaces each occurrence of SEARCH with REPLACE in FILE. Sed regular expressions are amazing but sometimes complicated so will be skipped here. [50 convenient commands](https://linuxhint.com/50_sed_command_examples/) or [more sed than you could ever want](https://www.grymoire.com/Unix/Sed.html)

### Generic one-liner for loops:

`for filepath in $(ls *fastq); do echo "$filepath" done`

Read this as "For each file that ends with fastq, store the path to that path as a variable named filepath, and print it"
* `$()` is an order-of-operations thing, anything inside `$()` is run first, in this case generating a list of files ending with "fastq"
* `for filepath in ` takes the list, and moves through the list, one at a time, and the currently read-in list entry is stored as the variable `filepath`
* `echo STRING` : print a string to the screen
* `$filepath` : the `$` means that the following word is a variable name, in this case the variable updated each iteration of the for loop

`for file in $(ls *fastq); do prefix=$(echo $file | sed 's/.fastq//'); mv $file $newname.fq; done`

* a simple renaming loop: for each file ending with "fastq", it defines a new variable called `prefix` by deleting the ".fastq" off with `sed` ("sample.fastq" would turned into "sample"), and then `mv` renames the original file to a new name, which is the contents of the prefix variable with a ".fq" suffix tacked on
* In bash, there can't be any spaces between variable name, equal sign, and contents. `prefix = 'something'` will be an error, but `prefix='SOMETHING` works

### Useful built-in variables
* `$HOME` : absolute path to your home folder
* `$PWD` : absolute path to the current directory
* `$PATH` : path(s) pointing to where shell should look for the program you want, like `grep`. If you're trying to download and add a tool you can run, make sure it's location is included in `$PATH`, either by moving it or updating `$PATH`

# Walkthrough / Tutorial Resources:
Some nice long-form resources:
[One](https://2017-dibsi-metagenomics.readthedocs.io/en/latest/command-line.html)    [Two](https://datacarpentry.org/shell-genomics/)    [Three](https://edcarp.github.io/shell-genomics-eddie/aio/)