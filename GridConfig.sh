#!/bin/bash
# File              : GridConfig.sh
# Author            : Anton Riedel <anton.riedel@tum.de>
# Date              : 25.08.2021
# Last Modified Date: 23.11.2021
# Last Modified By  : Anton Riedel <anton.riedel@tum.de>

# example configuration for running a analysis on grid

# check path to grid utility scripts which is ideally set in your bashrc or by sourcing SetupEnv.sh
# if not, set it to the default value and hope for the best
[ ! -d ${GRID_UTILITY_SCRIPTS:=$HOME/GridUtilityScripts} ] && echo "Cannot find GRID_UTILITY_SCRIPTS. Bail out..." && return 1

# miscellaneous variables
export TIMEOUT="300"
export COPY_JOBS="120"
export MAX_RESUBMIT="3"

# configure task
export TASK_BASENAME="Validation"
export ANALYSIS_NAME="Validation0"
export ALIPHYSICS_TAG="vAN-20210914_ROOT6-1"
export OUTPUT_TDIRECTORY_FILE="OutputAnalysis"
export INPUT_FILES_PER_SUBJOB="50"
export RUNS_PER_MASTERJOB="1"
export MASTER_RESUBMIT_THRESHOLD="50"
export TIME_TO_LIVE="44000"

# define directories and files on grid
export GRID_HOME_DIR="/alice/cern.ch/user/a/ariedel"
export GRID_WORKING_DIR_REL="${ALIPHYSICS_TAG}_${ANALYSIS_NAME}"
export GRID_WORKING_DIR_ABS="${GRID_HOME_DIR}/${GRID_WORKING_DIR_REL}"
export GRID_OUTPUT_DIR_REL="output" # relative to working directory on grid
export GRID_OUTPUT_DIR_ABS="${GRID_WORKING_DIR_ABS}/${GRID_OUTPUT_DIR_REL}"
export GRID_OUTPUT_ROOT_FILE="AnalysisResults.root"
# path to own xml collection so we do not read to recreate them every time
export GRID_XML_COLLECTION="${GRID_HOME_DIR}/XMLcollections/LHC10h/pass2/AOD/AOD160"
# analysis jdl file and analysis macro
export JDL_FILE_NAME="flowAnalysis.jdl"
export ANALYSIS_MACRO_FILE_NAME="flowAnalysis.C"

# define directories and files on local machine
export LOCAL_WORKING_DIR="$(realpath $(dirname ${BASH_SOURCE[0]}))"
export LOCAL_TMP_DIR="/tmp"
export LOCAL_OUTPUT_ROOT_FILE="Output.root"
export LOCAL_OUTPUT_HISTOGRAMS=$(
    cat <<'EOF'
[kRECO]fEventControlHistograms[kCEN][kBEFORE]
[kRECO]fEventControlHistograms[kCEN][kAFTER]
EOF
)

# set analysis mode
# has to be local or grid
export ANALYSIS_MODE="local"

# when runnnig on grid
# has to be test, offline, submit, full or terminate
# for submission, set runmode to offline
export GRID_RUN_MODE="offline"

# when running over real data, set to 1
# when running over monte carlo data, set to 0
export RUN_OVER_DATA="1"
# export RUN_OVER_DATA="0"
# when running over AOD, set to 1
# when running over ESD, set to 0
export RUN_OVER_AOD="1"
# export RUN_OVER_AOD="0"

[ $RUN_OVER_AOD -ne 1 ] && echo "Run over ESD not supported yet" && return 1

# when running locally specify data directory
[ $RUN_OVER_DATA -eq 1 -a $RUN_OVER_AOD -eq 1 ] && export DataDir="$HOME/data/2010/LHC10h/000137161/ESDs/pass2/AOD160"
[ $RUN_OVER_DATA -eq 0 -a $RUN_OVER_AOD -eq 1 ] && export DataDir="$HOME/sim/LHC10d4/120822/AOD056"

# configure centrality bins
export CENTRALITY_BIN_EDGES=$(
    cat <<'EOF'
0
5
EOF
)

# export CENTRALITY_BIN_EDGES=$(
#     cat <<'EOF'
# 0
# 5
# 10
# 20
# 30
# 40
# 50
# 60
# 70
# 80
# EOF
# )

# configure run numbers
export RUN_NUMBER=$(
    cat <<'EOF'
137161
EOF
)
# export RUN_NUMBER=$(
#cat <<'EOF'
# 139510
# 139507
# 139505
# 139503
# 139465
# 139438
# 139437
# 139360
# 139329
# 139328
# 139314
# 139310
# 139309
# 139173
# 139107
# 139105
# 139038
# 139037
# 139036
# 139029
# 139028
# 138872
# 138871
# 138870
# 138837
# 138732
# 138730
# 138666
# 138662
# 138653
# 138652
# 138638
# 138624
# 138621
# 138583
# 138582
# 138578
# 138534
# 138469
# 138442
# 138439
# 138438
# 138396
# 138364
# 138275
# 138225
# 138201
# 138197
# 138192
# 138190
# 137848
# 137844
# 137752
# 137751
# 137724
# 137722
# 137718
# 137704
# 137693
# 137692
# 137691
# 137686
# 137685
# 137639
# 137638
# 137608
# 137595
# 137549
# 137546
# 137544
# 137541
# 137539
# 137531
# 137530
# 137443
# 137441
# 137440
# 137439
# 137434
# 137432
# 137431
# 137430
# 137243
# 137236
# 137235
# 137232
# 137231
# 137230
# 137162
# 137161
# EOF
#)

echo "Sourced $(realpath ${BASH_SOURCE[0]}) at $(date "+%Y-%m-%d_%H:%M:%S")"

return 0
