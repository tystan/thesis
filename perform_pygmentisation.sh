



# to run type: 
# cd ~/Documents/thesis/R
# chmod 755 0000_getpygmentscode.sh
# ./0000_pygmentise_code.sh


cd ~/Documents/thesis/R

THISSTYLE='manni'
echo ${THISSTYLE}

pygmentize -f tex -S $THISSTYLE -a .syntax > tex/$THISSTYLE.tex


#######################################
## set variables 

TEX_OPTS0="mathescape=True,style="
THISSTYLE="autumn"
TEX_OPTS=$TEX_OPTS0$THISSTYLE
echo $TEX_OPTS
VERB_OPTS="verboptions=gobble=0,numbers=left,fontfamily=fvm,fontshape=n,fontsize=\footnotesize,tabsize=2"


#######################################
## use variables to set constants



#######################################
## these are the files to loop over

declare -a R_FILES=( \
'tophatcode' 'erosionquick' 'ctserosionquick' \
'quantnorm' 'MAadj' 'NMalign' 'guideTreePeakAlign' \
'doSVA' 'createPFldaobj' 'PFldapredict' \
'paretofronts' \
 );


#################################################
################## NEW: run this ################
#################################################

for t in "${R_FILES[@]}"
do
echo $t
pygmentize -O $TEX_OPTS -f tex -P $VERB_OPTS -o tex/$t".tex" R/$t".R"
done


pygmentize -O $TEX_OPTS -f tex -P $VERB_OPTS -o tex/domfeat.tex R/domfeat.c





