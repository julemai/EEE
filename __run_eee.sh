#!/bin/bash

# Copyright 2019 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of the EEE code library for "Computationally inexpensive identification
# of noninformative model parameters by sequential screening: Efficient Elementary Effects (EEE)".
#
# The EEE code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The MVA code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with The EEE code library.
# If not, see <https://github.com/julemai/EEE/blob/master/LICENSE>.
#
# If you use this method in a publication please cite:
#
#    M Cuntz & J Mai et al. (2015).
#    Computationally inexpensive identification of noninformative model parameters by sequential screening.
#    Water Resources Research, 51, 6417-6441.
#    https://doi.org/10.1002/2015WR016907.
#
set -e
#
# Perform a cleanup if script is interupted
trap cleanup 1 2 3 6
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$
#
# ---------------------------------------------------------------------------------------------------------------------
# ./run_eee.sh -p -c 1.24562587653 iter_1/
#
function usage () {
    printf "${pprog} [directory]                                                                                                \n"
    printf "Runs EEE in directory.                                                                                              \n"
    printf "                                                                                                                    \n"
    printf "Input                                                                                                               \n"
    printf "    directory        Directory containing file with masked parameters.                                              \n"
    printf "                                                                                                                    \n"
    printf "Options                                                                                                             \n"
    printf "    -h                    Prints this help screen.                                                                  \n"
    printf "    -m maskfile           Name of file containing information about model parameters (default: parameters.dat).     \n"
    printf "    -x model_function     Name of script that runs the model.                                                       \n"
    printf "                          (default: '2_run_model_ishigami-homma.py')                                                \n"
    printf "    -s modeloutputkey     Which model output will be analysed. Needs to be one of the keys used for                 \n"
    printf "                          dictionary of model outputs in '2_run_model_<name-model>.py'.                             \n"
    printf "                          (default: 'All').                                                                         \n"
    printf "                                                                                                                    \n"
    printf "Example                                                                                                             \n"
    printf "    ${isdir}/${pprog} -s out1 -x 2_run_model_ishigami-homma.py -m parameters.dat examples/ishigami-homma/           \n"
}
#
# cleanup at end wnd when interupted
function cleanup ()
{
  \rm -f *.${pid}
}
# ---------------------------------------------------------------------------------
# (0) Checking Arguments and Optionals
# ---------------------------------------------------------------------------------
#
# switches
maskfile='parameters.dat'
outfile='eee_results.pdf'
cutoff=-1
traj_M1=5    # number of trajectories for first iteration
traj_M2=5    # number of trajectories for final iteration
traj_M=1     # number of trajectories for 2nd, 3rd, ..., second-last iteration
model_function='2_run_model_ishigami-homma.py'
modeloutputkey='All'

verbose=2 # 0: pipe stdout and stderr to /dev/null
          # 1: pipe stdout to /dev/null
          # 2: print everything
pipeit=''
if [[ ${verbose} -eq 0 ]] ; then pipeit=' > /dev/null 2>&1' ; fi
if [[ ${verbose} -eq 1 ]] ; then pipeit=' > /dev/null' ; fi

while getopts "hpm:s:x:" Option ; do
    case ${Option} in
        h) usage 1>&2; exit 0;;
        m) maskfile="${OPTARG}";;
        s) modeloutputkey="${OPTARG}";;
        x) model_function="${OPTARG}";;
        *) printf "Error ${pprog}: unimplemented option.\n\n";  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))
#
# Check args
NO_ARGS=1
if [[ $# -lt ${NO_ARGS} ]] ; then
    printf "Error ${pprog}: directory missing.\n\n"
    usage 1>&2
    exit 1
fi
indir="$@"
#
# Check input directory exist
if [ ! -d ${indir} ] ; then
    printf "Error ${pprog}: Input directory not found: %s\n" ${indir}
    exit 1
fi
#
# Check if model output key is given
if [[ ${modeloutputkey} == 'None' ]] ; then
    printf "Error ${pprog}: Model output key must be given. Needs to be one of the keys used for model output dictionary in 2_run_model_<name-model>.py. \n"
    exit 1
fi

#
# Change to given directory
cd ${indir}
#
# Check if all necessary files are in directory
if [ ! -f ${maskfile} ] ; then
    printf "Error ${pprog}: Parameter information file not found: %s/%s\n" ${indir} ${maskfile}
    exit 1
fi

first_iteration=true    # true for first iteration
last_iteration=false    # true if last iteration is running (after two consecutive iterations lead to same number of non-informative parameters)
finished=false          # true when last iteration is finished
iterations_counter=1    # number of iterations of EEE
n_model_runs=0          # number of model runs performed in total
n_informative=0         # number of informative parameters

# ------------------------------------------
# Here starts the loop for iterations
# ------------------------------------------
while [[ "${finished}" = false ]] ; do

    # Create new folder for next iteration
    mkdir iter_${iterations_counter}

    # Get mask_file (either initial one           OR from iteration before) and
    #     cut-off   (either -1=not determined yet OR from cutoff.dat file from iteration before)
    if [ ${iterations_counter} -eq 1 ] ; then
        cp ${maskfile} iter_${iterations_counter}/.
        cutoff=-1
    else
        iter_before=$[${iterations_counter}-1]
        cp iter_${iter_before}/${maskfile}.new iter_${iterations_counter}/${maskfile}
    fi

    # Determine number of trajectories
    if ${first_iteration} ; then
        # How many trajectories at beginning?
        traj=${traj_M1}
    else
        if ${last_iteration} ; then
        # How many trajectories at end?
        traj=${traj_M2}
        else
        # How many trajectories in between?
        traj=${traj_M}
        fi
    fi

    # Change directory to iteration-folder
    cd iter_${iterations_counter}

    echo '# ---------------------------------------------------------------------------------'
    echo '# ('${iterations_counter}'.1) Create Morris trajectories                           '
    echo '# ---------------------------------------------------------------------------------'
    python "${isdir}"/codes/1_create_parameter_sets.py -d "${maskfile}" -t ${traj} -n 1 -o parameter_sets

    echo '# ---------------------------------------------------------------------------------'
    echo '# ('${iterations_counter}'.2) Run model and store all model results                  '
    echo '# ---------------------------------------------------------------------------------'
    parafile_M=$( \ls parameter_sets_1_scaled_*_M.dat )
    parafile_v=$( \ls parameter_sets_1_*_v.dat )
    skip=$( head -1 ${parafile_M} | cut -d : -f 2 )                # number of lines to skip in file containing sampled parameter sets
    nlines=$( echo $( wc -l ${parafile_M}) | cut -f 1 -d " ")      # total number of lines in file containing sampled parameter sets
    n_model_runs=$(( ${n_model_runs} + ${nlines} - ${skip} ))      # number of model runs
    echo 'number model runs: '${n_model_runs}

    python "${isdir}"/codes/${model_function} -i "${parafile_M}" -s ${skip} -o model_output.pkl

    echo '# ---------------------------------------------------------------------------------'
    echo '# ('${iterations_counter}'.3) Calculate Elementary Effects                         '
    echo '# ---------------------------------------------------------------------------------'
    model_outputs='model_output.pkl'
    eefile='eee_results.dat'
    parafile_M=$( \ls parameter_sets_1_*_M.dat | grep -v scaled )
    parafile_v=$( \ls parameter_sets_1_*_v.dat )
    python "${isdir}"/codes/3_derive_elementary_effects.py -i ${model_outputs} -k ${modeloutputkey} -d "${maskfile}" -m "${parafile_M}" -v "${parafile_v}"  -o "${eefile}"

    echo '# ---------------------------------------------------------------------------------'
    echo '# ('${iterations_counter}'.4) Create some plots and derive cutoff                  '
    echo '# ---------------------------------------------------------------------------------'
    python "${isdir}"/codes/4_derive_threshold.py -e "${eefile}" -m "${maskfile}" -p "${outfile}" -c ${cutoff} #-t # -n

    echo '# ---------------------------------------------------------------------------------'
    echo '# ('${iterations_counter}'.5) Get the Cutoff                                       '
    echo '# ---------------------------------------------------------------------------------'
    files=$(\ls cutoff_*.dat)
    cutoff=''
    for ff in $files ; do
        cutoff=$(echo ${cutoff}$(tail -1 $ff):)
    done
    echo 'Cutoff(s) :: '${cutoff}
    if ${first_iteration} ; then
        first_iteration=false
    fi

    #
    # Grep flags:
    # -F  Interpret PATTERN as a list of fixed strings, separated by newlines, any of which is to be matched.
    # -x  Select only those matches that exactly match the whole line.
    # -v  Invert the sense of matching, to select non-matching lines.
    # -f  Obtain patterns from FILE, one per line.  The empty file contains zero patterns, and therefore matches nothing.
    # new_parameters_detected=$(grep -Fxvf "iter_${iterations_counter}/mask_para.dat" "iter_${iterations_counter}/mask_para.dat.new" | wc -l)
    # echo '      In this iteration '${new_parameters_detected}' parameters where additionally detected to be informative'
    # n_informative=$[${n_informative}+${new_parameters_detected}]

    # Compare 2nd columns (transformation_type) in mask_para.dat and mask_para.dat.new
    # to determine
    # (1) number of additionally detected informative parameters :: new_parameters_detected
    # (2) number of non-informative parameters left              :: n_non_informative
    new_parameters_detected=0
    lines=$(echo $(wc -l ${maskfile}.new) | cut -f 1 -d " ") # new files has no blank lines at the end
    n_non_informative=0
    for ((ii=3 ; ii<=${lines} ; ii++)) ; do

        mask_old=$(echo $(head -${ii} ${maskfile}     | tail -1) | cut -f 6 -d " ")
        mask_new=$(echo $(head -${ii} ${maskfile}.new | tail -1) | cut -f 6 -d " ")

        if (( ${mask_old} != ${mask_new} )) ; then
            new_parameters_detected=$((new_parameters_detected+1))
        fi
        if (( ${mask_new} > 0 )) ; then
            n_non_informative=$((n_non_informative+1))
        fi
    done
    echo 'In this iteration '${new_parameters_detected}' parameters where additionally detected to be informative.'
	echo "n informative before = ${n_informative}"
    n_informative=$[${n_informative}+${new_parameters_detected}]

    if ${last_iteration} ; then
        finished=true
    fi

    # stop here since no non-informative parameter is left
    if (( ${n_non_informative} == 0 )) ; then
        finished=true
    fi

    if [[ "${finished}" = false ]] ; then
        # Determine if new and old mask_file are the same
        # Yes --> Next iteration is the last one
        # No  --> Go one with an intermediate iteration
        if [ ${new_parameters_detected} -gt 0 ] ; then
            last_iteration=false
        else
            last_iteration=true
        fi
        iterations_counter=$[${iterations_counter}+1]
    fi

    cd ..

	echo "New parameters detected =  ${new_parameters_detected}"
	echo "n informative after = ${n_informative}"

done # Loop over iterations as long as ${finished} is not true

cp "iter_${iterations_counter}/${maskfile}.new" ${maskfile}.final

touch "eee_info.dat"
echo "Info about EEE analysis"                                   >  "eee_info.dat"
echo "------------------------------------------------"          >> "eee_info.dat"
echo "Number of informative parameters :: "${n_informative}      >> "eee_info.dat"
echo "Number of model runs in total    :: "${n_model_runs}       >> "eee_info.dat"
echo "Number of iterations             :: "${iterations_counter} >> "eee_info.dat"

echo ''
echo '# ---------------------------------------------------------------------------------'
echo "# Info about EEE analysis"
echo '# ---------------------------------------------------------------------------------'
echo "      Number of informative parameters :: "${n_informative}
echo "      Number of model runs in total    :: "${n_model_runs}
echo "      Number of iterations             :: "${iterations_counter}
echo "      EEE info file written to "${indir}"eee_info.dat"
echo " "
echo " "

exit 0
