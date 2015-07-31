



### to run at once in the bash shell: 
# chmod 755 perform_pygmentisation.sh
# ./perform_pygmentisation.sh

# set curr dir
cd ~/Documents/ThesisAppRcode

#######################################
## set variables 


THISSTYLE="manni" # doesn't matter for generating 
TEX_OPTS="mathescape=True,style="$THISSTYLE
echo $TEX_OPTS
VERB_OPTS="verboptions=gobble=0,numbers=left,fontfamily=fvm,fontshape=n,fontsize=\footnotesize,tabsize=2"

#######################################
## style file:
## uncomment to produce the style file

# THISSTYLE='manni'
# echo ${THISSTYLE}
# pygmentize -f tex -S $THISSTYLE -a .syntax > tex/$THISSTYLE.tex

#######################################
## these are the files to loop over

# create array of R script filenames
cd R
declare -a R_FILES=(*.R)
cd ..
echo ${R_FILES[@]}
# length of array
J=${#R_FILES[@]}
echo $J

for ((j=0; j<$J; j++))
do
    #remove the ".R" from each string/filename
    R_FILES[$j]=${R_FILES[$j]%.R}
    echo "${R_FILES[$j]}"
done


#############################################
## use pygments on .R files --> .tex files 

for t in "${R_FILES[@]}"
do
echo $t
pygmentize -O $TEX_OPTS -f tex -P $VERB_OPTS -o tex/$t".tex" R/$t".R"
done

# also the one C file
pygmentize -O $TEX_OPTS -f tex -P $VERB_OPTS -o tex/11_dom_feat.tex R/11_dom_feat.c





