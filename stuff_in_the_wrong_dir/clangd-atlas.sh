#!/bin/bash

# Build the full clangd command string with all arguments quoted properly
quoted_args=""
for arg in "$@"; do
    quoted_args+=" \"$arg\""
done


# Build the postsetup command
cmd="/home/padi617e/.vscode-server/data/User/globalStorage/llvm-vs-code-extensions.vscode-clangd/install/19.1.2/clangd_19.1.2/bin/clangd${quoted_args}"

# Source the ATLAS environment and run clangd in postsetup
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase;
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh -c el9 --runpayload="$cmd"
