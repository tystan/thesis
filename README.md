## PhD dissertation appendix files 
@ [`github.com/tystan/thesis/`](https://github.com/tystan/thesis/)

R functions and code for the PhD dissertation of Tyman Stanford.

> Statistical analysis of proteomic mass spectrometry data for the identification of biomarkers and disease diagnosis. The University of Adelaide, 2015.


## Navigation

* A PDF version of the code is available in [`AppendixRCode.pdf`](../master/AppendixRCode.pdf). 
    * This contains all the relevant R code in one syntax highlighted self-contained document.
    * [`AppendixRCode.tex`](../master/AppendixRCode.tex) top-level of this repository is the TeX file required to generate [`AppendixRCode.pdf`](../master/AppendixRCode.pdf). 
    * The directory [`tex/`](../../tree/master/tex/) contains the source files referenced in [`AppendixRCode.tex`](../master/AppendixRCode.tex). 
* The individual R files are available in the [`R/`](../../tree/master/R/) directory. 
    * They are numbered so as to be ordered in sequence of appearance in the dissertation and [`AppendixRCode.pdf`](../blob/master/AppendixRCode.pdf).


## Acknowledgments

The wonderful [Pygments](http://pygments.org/) command-line tool was used to peform syntax highlighting. [Pygments](http://pygments.org/) can convert `R` files (and many other programming languages) into `Verbatim` environment `TeX` files to display syntax highlighted code in `LaTeX` documents. Many thanks to this freely available software.
* The file [`perform_pygmentisation.sh`](../master/perform_pygmentisation.sh) is the command-line input (for unix-alike platforms) used to convert the files in [`R/`](../../tree/master/R/) to the files in [`tex/`](../../tree/master/tex/).

##### Installing pygments

I was able to install pygments with the following steps: 
```sh
#### download `get-pip.py' from: https://bootstrap.pypa.io/get-pip.py

#### then in terminal/shell, move to the location of get-pip.py
cd ~/{dir-location-of-get-pip.py}
#### and run get-pip.py in python: (use `sudo' for admin access)
sudo python get-pip.py

#### now to install pygments, type:
sudo pip install pygments
#### how'd it go? (there should be files listed here)
pip show --files pygments
```
A more thourough treatment is available at the download pages of [pygments](http://pygments.org/download/) and [pip](https://pip.pypa.io/en/stable/installing.html) should you run into any trouble.
