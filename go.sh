source activate py2.7

Folder=/projectZ/tseng/milkyway/data

LAST=`ls -l ${Folder} | wc -l`
LAST=$(printf "%02d" $LAST)

mkdir                    $Folder/$LAST


python dump_file.py |& tee $Folder/$LAST/log
conda deactivate

mv UM_IC                 $Folder/$LAST
mv ExtPotTable           $Folder/$LAST
mv FractalDensity       $Folder/$LAST
mv FractalUxyz          $Folder/$LAST
cp *.py                  $Folder/$LAST 
echo "Done"
