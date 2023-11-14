#!/bin/bash
helpFunction()
{
   echo "verabole i is the in pdb file"
   echo "verabole or is the organisem to sharch aginst"
   echo "verabole o is the out file name"
   echo "for example bash usalign_af.sh -i /home/nirc/usalign/pdb_files/OG1032.txt -a all -o OG1032 "
   echo "for example bash usalign_af.sh -a group -o try1 -c parallel -in yes"
   exit 1 # Exit script after printing help
} 
#the part above is not relevant ^
#cluster & scedualer type(maybe as configfile?) or local run
#out type add msa and stracture suprepos

while getopts "i:a:o:c:sn:in:" opt
do
   case "$opt" in
      i ) in_file="$OPTARG" ;;
      a ) batch_type="$OPTARG" ;;
      o ) out_syntax="$OPTARG" ;;
      c ) cluster="$OPTARG" ;;
      sn ) standart_AF_names="$OPTARG" ;;
      in ) indexing="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done
if [ -z "$batch_type" ] || [ -z "$out_syntax" ]; then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

## verabole to add for other use
exet = "USalign"
numbee_of_cores_limit  = 16



if [ -z "$cluster" ]; then
   cluster="ST"
fi

# if [ $cluster == "bsub" ]; then
#    log_file_location=$home_folder/log/$out_syntax/
#    error_file_location=$home_folder/err/$out_syntax/
#    mkdir -p $log_file_location
#    mkdir -p $error_file_location
#    scedualer_submition_commend="bsub -q new-long -o $log_file_location -e $error_file_location"
# elif [ $cluster == "slurm" ]; then
# # add info
#    log_file_location=$home_folder/log/$out_syntax/
#    error_file_location=$home_folder/err/$out_syntax/
#    mkdir -p $log_file_location
#    mkdir -p $error_file_location
#    scedualer_submition_commend="bsub -q new-long -o $log_file_location -e $error_file_location"
elif [ $cluster == "parallel" ]; then
   CHECK_MARK="\033[0;32m\xE2\x9C\x94\033[0m"
   counter=1
fi

home_folder=/home/nirc/usalign
stracture_folder=/home/nirc/usalign/pdb_files/

if [[ "$indexing" = "yes" ]]; then
#|| [ $indexing == "y" ]; then
   foldlist=`find $stracture_folder -maxdepth 1 -type d`
   for folder in $foldlist; do
      orto_name=${folder##*/}
      ls $folder/*.pdb | xargs -n1 -I{} basename "{}">$stracture_folder$orto_name.txt
   done
   #statements
fi
## tis part process the argument for the submition 
## in this case its based on the list text files with paths for bdp files
# so you will need to re make this part 
if [ $batch_type == "all" ]; then
   batch_list=()
   list_of_stractures_list=`ls $stracture_folder*.txt`
   for list_of_stractures in $list_of_stractures_list; do
      batch_list+=(${list_of_stractures##*/})
   done
   #extract protein id
   #remove full path from file name
   filename=${in_file##*/}
   protein_id=${filename%.*}
   # makeing the folder 
   mkdir -p $home_folder/$out_syntax
   # building the array stracures
   mapfile -t in_file_map < $in_file
   # this option was added to deal with AF predicted protein naming you will need to replace it
   for p2_stracture_file_full_path in ${in_file_map[*]}; do
      if [ $standart_AF_names == "YES" ] || [ $standart_AF_names == "y" ]; then   
         p2_stracture_file=${p2_stracture_file_full_path##*/}
         p2_file_name_with_ext=${p2_stracture_file%.*}
         p2_id_with_AF_id=${p2_file_name_with_ext%-F*-model_??}
         p2_prot_id_array=${p2_id_with_AF_id##"AF-"}
         p2_prot_id=($p2_prot_id_array)
      else
         p2_stracture_file=${p2_stracture_file_full_path##*/}
         p2_prot_id_array=${p2_stracture_file%.*}
         p2_prot_id=($p2_prot_id_array)
      fi
      elif [ $cluster == "parallel" ]; then
         for file_in_list in ${batch_list[*]}; do
            folder_to_run_against=${file_in_list%.*.*}
            if [[ "$folder_to_run_against" == *"."* ]]; then
               folder_to_run_against=${file_in_list%.*}
            fi
            # this is the submition part 
            # $home_folder/$exet this is the exetacutabole folowing the arguments
            # the importent part is the & at the end to submit as a sub process 

            $home_folder/$exet -split 0 $stracture_folder$protein_id/$p2_stracture_file_full_path -dir2 $stracture_folder$folder_to_run_against/ $stracture_folder$file_in_list -outfmt 2 > $home_folder/$out_syntax/$p2_prot_id_1-$file_in_list&
            echo "runing orto group $p2_prot_id_1"
         done 
         # check number of processes with the exet name
         if [ $(pgrep -c $exet) -lt $numbee_of_cores_limit ]; then
            t=1
         else
            # if its more then the limit its got to a while loop until its eq
            while  [ $(pgrep -c $exet) -eq $numbee_of_cores_limit ]; do
               echo -n "$(pgrep -c $exet) runing processes, waiting!"
               sleep 5s
               echo -e -n "\\r\033[K"
            done
            echo -e  "\\r${CHECK_MARK} $counter  processes done"
            let counter++
         fi 