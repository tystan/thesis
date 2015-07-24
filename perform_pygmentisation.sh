



### to run at once: 
# chmod 755 0000_getpygmentscode.sh
# ./0000_pygmentise_code.sh


cd ~/Documents/ThesisAppRcode

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
'00_erosion_slow' '05_ma_adj' '10_pfda_predict' \
'01_erosion_quick' '06_create_w' \
'02_cts_erosion_slow' '07_dendro_peak_align' '12_pareto_fronts' \
'03_cts_erosion_quick' '08_do_sva' \
'04_quant_norm' '09_create_pfda_obj' \
 );
echo ${R_FILES[@]}

#################################################
################## NEW: run this ################
#################################################

for t in "${R_FILES[@]}"
do
echo $t
pygmentize -O $TEX_OPTS -f tex -P $VERB_OPTS -o tex/$t".tex" R/$t".R"
done


pygmentize -O $TEX_OPTS -f tex -P $VERB_OPTS -o tex/11_dom_feat.tex R/11_dom_feat.c





