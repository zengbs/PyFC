source activate py2.7

Folder=/projectZ/tseng/milkyway/data

LAST=`ls $Folder|sort -n|tail -n1`
LAST=`python -c "print(int(\"$LAST\")+1)"`
LAST=$(printf "%03d" $LAST)

mkdir                      $Folder/$LAST

python dump_file.py |& tee $Folder/$LAST/log
conda deactivate

python double2float.py

mv UM_IC*                  $Folder/$LAST
mv ExtPotTable*            $Folder/$LAST
mv FractalDensity*         $Folder/$LAST
mv FractalUx*              $Folder/$LAST
mv FractalUy*              $Folder/$LAST
mv FractalUz*              $Folder/$LAST
cp *.py                    $Folder/$LAST 
echo "Done"
