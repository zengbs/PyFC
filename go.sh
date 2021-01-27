
LAST=`ls -l /projectZ/tseng/milkyway/data/ | wc -l`
LAST=$(printf "%02d" $LAST)

mkdir /projectZ/tseng/milkyway/data/$LAST
source activate py2.7
python dump_file.py |& tee /projectZ/tseng/milkyway/data/$LAST/log
conda deactivate

mv UM_IC                 /projectZ/tseng/milkyway/data/$LAST
mv ExtPotTable           /projectZ/tseng/milkyway/data/$LAST
mv FractalDensity        /projectZ/tseng/milkyway/data/$LAST
mv FractalUxyz           /projectZ/tseng/milkyway/data/$LAST
cp *.py                  /projectZ/tseng/milkyway/data/$LAST 
echo "Done"
