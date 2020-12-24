#!/bin/bash

rm -rf UM_IC ExtPotTable
mv parameters.py parameters_temp.py
cp regression_test/parameters_example.py parameters.py
python dump_file.py
sha1sum Mukherjee_UM_IC       regression_test/Standard_UM_IC
sha1sum Mukherjee_ExtPotTable regression_test/Standard_ExtPotTable
rm parameters.py
mv parameters_temp.py parameters.py

