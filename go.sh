source activate py2.7
python dump_file.py
conda deactivate

LAST=`ls -l data/ | wc -l`
LAST=$(printf "%02d" $LAST)

mkdir data/$LAST
mv Mukherjee_UM_IC       data/$LAST
mv Mukherjee_ExtPotTable data/$LAST
cp parameters.py         data/$LAST 
echo "Done"
