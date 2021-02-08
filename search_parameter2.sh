#!/bin/bash

GetJobID()
{
  local stdout=`$@` # call the command and capture the stdout
  local id=`printf  "$stdout" | awk -F. '{print $1}'` # get the jobid
  echo  "$id"
}

WaitForHULK()
{
   shopt -s expand_aliases

   local sleep_time=1 # seconds; don't make this too short! don't want to tax system with excessive qstat calls
   local me=`whoami`
   local ID
   alias myqstat='qstat | grep $me'

   for ID in $@
   do
      local S=`myqstat | grep $ID | awk '{print $5}'` # check to see if job is running
	         while [[ "$S" == "R" || "$S" == "Q" ]] # while $status is runing or in qune
		 do
			 sleep $sleep_time
			 S=`myqstat | grep $ID | awk '{print $5}'`
		 done
       echo "Job ID:$ID is done!"
   done
}

Uxyz_kmin=1.0
Uxyz_mean=1.0
Uxyz_sigma=1e6


sed -i  "s/^  Uxyz_kmin.*/  Uxyz_kmin     = ${Uxyz_kmin}/"                parameters.py
sed -i  "s/^  Uxyz_mean.*/  Uxyz_mean     = ${Uxyz_mean} \/ Const_C/"     parameters.py
sed -i  "s/^  Uxyz_sigma.*/  Uxyz_sigma    = ${Uxyz_sigma} \/ Const_C/"   parameters.py


qsub submit.job


JobID=`ls JobID.*`
ID="${JobID##*.}"


FOLDER=Uxyz_kmin_${Uxyz_kmin}_Uxyz_mean_${Uxyz_mean}_Uxyz_sigma_${Uxyz_sigma}

mkdir /projectZ/tseng/jet-turning/Yang-frac_dens-o-frac_vel-o-50kpc-gravity-o-jet-x/${FOLDER}
cd    /projectZ/tseng/jet-turning/Yang-frac_dens-o-frac_vel-o-50kpc-gravity-o-jet-x/${FOLDER}

rsync -avzhP /projectZ/tseng/jet-turning/Yang-frac_dens-o-frac_vel-o-50kpc-gravity-o-jet-x/temp/* .

LAST=`ls /projectZ/tseng/milkyway/data/ |sort -n|tail -n1`

ln -fs /projectZ/tseng/milkyway/data/${LAST}/UM_IC .
ln -fs /projectZ/tseng/milkyway/data/${LAST}/ExtPotTable .


cd /projectZ/tseng/jet-turning/Yang-frac_dens-o-frac_vel-o-50kpc-gravity-o-jet-x/${FOLDER}

WaitForHULK $ID
qsub submit.job
