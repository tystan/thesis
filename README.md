## :blue_book: PhD dissertation online reference: selected appendices containing `R` code   :blue_book:
@ [`github.com/tystan/thesis/`](https://github.com/tystan/thesis/)

`R` functions and code for the PhD dissertation of Tyman Stanford.

> Statistical analysis of proteomic mass spectrometry data for the identification of biomarkers and disease diagnosis. The University of Adelaide, 2015.


## :orange_book: Navigation :orange_book:

* A PDF version of the `R` code used and referenced throughout the dissertation is available in [`AppendixRCode.pdf`](../master/AppendixRCode.pdf). 
    * This contains all the relevant `R` code in one syntax highlighted self-contained document.
* The individual `R` files are available in the [`R/`](../../tree/master/R/) directory. 
    * The files are numbered so as to be ordered in sequence of appearance in the dissertation and [`AppendixRCode.pdf`](../blob/master/AppendixRCode.pdf).
* The top-level file [`AppendixRCode.tex`](../master/AppendixRCode.tex) is the TeX file required to generate [`AppendixRCode.pdf`](../master/AppendixRCode.pdf). 
    * The directory [`tex/`](../../tree/master/tex/) contains the source files referenced in [`AppendixRCode.tex`](../master/AppendixRCode.tex). 


## :closed_book: Acknowledgments :closed_book:

The [Pygments](http://pygments.org/) command-line tool was used to perform syntax highlighting. [Pygments](http://pygments.org/) can convert `R` files (and many other programming languages) into syntax highlighted `Verbatim`environment `LaTeX` markup as `.tex` files. Many thanks to this freely available software.
* The file [`perform_pygmentisation.sh`](../master/perform_pygmentisation.sh) is the command-line input (for unix-alike platforms) used to convert the files in [`R/`](../../tree/master/R/) to the files in [`tex/`](../../tree/master/tex/). This requires Pygments to be installed.

##### Installing pygments

The following steps were required to install [Pygments](http://pygments.org/): 
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
A more thorough treatment is available at the download pages of [Pygments](http://pygments.org/download/) and [pip](https://pip.pypa.io/en/stable/installing.html) should it be required.
