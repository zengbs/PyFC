source activate py2.7

Folder=/projectZ/tseng/milkyway/data

LAST=`ls -l ${Folder} | wc -l`
LAST=$(printf "%02d" $LAST)



python dump_file.py |& tee $Folder/$LAST/log
#python dump_file.py
conda deactivate

mkdir                    $Folder/$LAST
mv UM_IC                 $Folder/$LAST
#mv ExtPotTable           $Folder/$LAST

if [ -f "FractalDensity" ]; then
   mv FractalDensity       $Folder/$LAST
fi

if [ -f "FractalUxyz" ]; then
   mv FractalUxyz          $Folder/$LAST
fi

cp *.py                  $Folder/$LAST 
echo "Done"
