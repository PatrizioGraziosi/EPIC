#!/bin/bash
#SBATCH --time=01:59:00
#SBATCH --job-name=DPsHH24

#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks-per-node=1    # Number of MPI ranks per node
#SBATCH --cpus-per-task=1      # Number of OpenMP threads for each MPI process/rank
#SBATCH --mem=7000           # Per nodes memory request (MB)
#SBATCH --account=IscrB_IPLAS-X

#SBATCH --partition=g100_usr_prod


n_q=6
n_modes=72
n_atoms=24
modes_to_skip=0
fileName_bxsf="'HH24.bxsf'"
material_name="'HH24'"
fileName="modes_info.m"
E_band_window=2 # in eV, range beyond the band edge where the bands are considered, it is related to what you need in transport simulations, for semicondcuctors it is easilt 0.5, for metals likely needs to be wider
labels_q="{'G', 'X', 'T', 'W', 'S', 'R' }" # q-point labels for the tables


total_modes=$((n_atoms*3))

printf " n_q = $n_q; \n n_modes = $n_modes; \n n_atoms = $n_atoms; \n total_modes = $total_modes; \n modes_to_skip = $modes_to_skip; \n fileName = $fileName_bxsf; \n material_name = $material_name; \n E_band_window = $E_band_window ; \n labels_q = $labels_q ; \n" > $fileName


module load profile/eng
# module load intel/oneapi-2021--binary
# module load intelmpi/oneapi-2021--binary
module load matlab/r2022b



sed -n "3,5p" ../MPOSCAR-orig > A_matrix.txt

line_POSCAR=$(( $n_atoms+8 ))
sed -n "9,${line_POSCAR}p" ../MPOSCAR-orig > coord_initial.txt


# the matlab script that generates the line for the setting....conf file runs, for each q-point, all the modes
counter=0
for iq in `seq 1 $n_q ` ; do
	for im in `seq 1 $n_modes ` ; do

        	counter=$(( $counter+1 ))
		if (( $counter <= 9 )); then
			mod_dir="mod-00$counter"
		elif (( $counter <= 99 )); then
                        mod_dir="mod-0$counter"
		else
			mod_dir="mod-$counter"
                fi
		cd "$mod_dir"
			cp ../bxsf_to_ELECTRA_shifting_option.m ./
			cp ../shifting_bands.m ./
                	sed -n "9,${line_POSCAR}p" POSCAR > coord.txt
                cd ../
	done
done
cd mod-all
     cp ../bxsf_to_ELECTRA_shifting_option.m ./
     cp ../shifting_bands.m ./
     sed -n "9,${line_POSCAR}p" POSCAR > coord.txt
cd ../


cp mod-orig/*bxsf ./

matlab -nosplash <bands_mod_extract_DefPot_narrowbands_shifted_script_L.m> DefPot.log
matlab -nosplash <average_DPs.m>DP.out
