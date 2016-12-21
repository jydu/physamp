'''!PhySamp''' is a package dedicated to phylogenetic sampling. It samples a sequence alignment according to its corresponding phylogenetic tree. Current version is 0.2.0.


# Description

The !PhySamp package currently contains one program, '''bppalnoptim''', which samples a sequence alignment by removing sequences in order to maximize the number of sites suitable for a given analysis.
The program has three running modes:
* Interactive: the user will be iteratively proposed a set of choices for sequence removal, with their corresponding site gains. The procedure stops when the user does not want to remove more sequences, and the resulting filtered alignment is written.
* Automatic: the user enters an a priori criterion for stopping the filtering procedure (for instance a minimum number of sequences to keep).
* Diagnostic: this mode allows to plot the trade-off curve, by showing the site gain as a function of the number of removed sequences.


# Availability

The '''bppalnoptim''' program is command-line driven. You can get pre-compiled executable files for your system (if there are any), use pre-compiled packages (if there are any) or compile the programs yourself (should work on any system with a decent C++ compiler).
The latest version of !PhySamp (0.2.0) is based on Bio++ 2.2.0 [[http://biopp.univ-montp2.fr/]].

 * Pre-compiled executables are available here [linux x64](http://biopp.univ-montp2.fr/repos/exe/lin64/physamp) [linux i386](http://biopp.univ-montp2.fr/repos/exe/lin32/physamp)
 * Source code is available here [here](http://biopp.univ-montp2.fr/repos/sources/physamp/)
 * Git repository address is https://github.com/jydu/physamp

The programs depend on the Bio++ libraries. Pre-compiled executables are statically linked, and therefore already include all required code from the libraries. Pre-compiled packages will ask for all required dependencies, which can be found in the same download directory. For compiling the programs yourself, from the downloaded sources or from the git repository, please follow the instructions from the Bio++ website http://biopp.univ-montp2.fr/wiki/index.php/Installation.

# Usage

Several example data sets are distributed along with the source code of the package. A reference manual is also available [here](http://biopp.univ-montp2.fr/manual/html/physamp/), or can be downloaded as [PDF](http://biopp.univ-montp2.fr/manual/pdf/physamp). Questions can be asked on the dedicated forum: [here](https://groups.google.com/forum/?hl=en#!forum/physamp).

