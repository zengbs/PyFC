source activate py2.7
python dump_file.py |& tee data/$LAST/log
conda deactivate

LAST=`ls -l data/ | wc -l`
LAST=$(printf "%02d" $LAST)

mkdir data/$LAST
mv UM_IC                 data/$LAST
mv ExtPotTable           data/$LAST
mv FractalDensity        data/$LAST
mv FractalUxyz           data/$LAST
cp *.py                  data/$LAST 
echo "Done"
