#!/bin/bash

#######################################################################
# Usage: submit.sh [--project=PROJECT] [--tag=TAG]
#
# Arguments:
#   --project=PROJECT   : Specify the project directory
#   --tag=TAG           : Git ref to checkout on grid nodes (default: develop)
#######################################################################

# Print usage information
usage() {
  echo "Usage: submit.sh [--project=PROJECT] [--tag=TAG]"
  echo ""
  echo "Arguments:"
  echo "  --project=PROJECT   : Specify the project directory"
  echo "  --tag=TAG           : Git ref to checkout on grid nodes (default: develop)"
}

# Initialize variables
PROJECT=""
TAG="ccpionproton"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project=*)
      PROJECT="${1#*=}"
      shift
      ;;
    --project)
      PROJECT="$2"
      shift 2
      ;;
    --tag=*)
      TAG="${1#*=}"
      shift
      ;;
    --tag)
      TAG="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --) # end of options
      shift
      break
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

#######################################################################
# Check for required arguments
#######################################################################
missing_args=()
[[ -z "$PROJECT" ]] && missing_args+=("--project")

if [[ ${#missing_args[@]} -gt 0 ]]; then
    echo "Error: Missing required argument(s): ${missing_args[*]}" >&2
    usage
    exit 1
fi

#######################################################################
# Initial Setup
#######################################################################

# IFDH options
export IFDH_CP_MAXRETRIES=0
export IFDH_WEB_TIMEOUT=100

COPY_ATTEMPTS=3
COPY_RETRY_DELAY_SECONDS=10

copy_to_local_with_retry() {
    local source_path="$1"
    local destination_path="$2"
    local description="$3"
    local attempt

    for ((attempt = 1; attempt <= COPY_ATTEMPTS; attempt++)); do
        # A failed gfal-copy can leave a zero-byte local destination. Remove
        # only this job-local file before retrying so that it cannot be
        # mistaken for a successful copy.
        rm -f -- "$destination_path"

        echo "Copying ${description} (attempt ${attempt}/${COPY_ATTEMPTS}): ${source_path}"
        if ifdh cp "$source_path" "$destination_path" && [[ -s "$destination_path" ]]; then
            return 0
        fi

        echo "Warning: failed to copy ${description} on attempt ${attempt}." >&2
        if ((attempt < COPY_ATTEMPTS)); then
            sleep "$COPY_RETRY_DELAY_SECONDS"
        fi
    done

    echo "Error: failed to copy ${description} after ${COPY_ATTEMPTS} attempts." >&2
    return 1
}



# Setup CVMFS area
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

# Setup the required dependencies
setup sbnana v10_01_02_01 -q e26:prof
setup cmake v3_27_4

ups active

# Build medulla
git clone https://github.com/kyjung123/medulla.git
cd medulla
git checkout ${TAG}
mkdir build && cd build
export CC=$(which gcc)
export CXX=$(which g++)
cmake .. -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC
make -j4

#######################################################################
# Prestage job-specific files
#######################################################################

# Copy the project database
#ifdh cp $PROJECT/project.db project.db
if ! copy_to_local_with_retry "$PROJECT/project.db" "project.db" "project database"; then
    exit 1
fi


# Extract this job's configuration file. First, we get the job ID for this
# process by checking against the list of not-yet-completed jobs in the project
# database.
JOBID=$(sqlite3 -noheader project.db "SELECT jobid FROM jobs WHERE status != 'completed' ORDER BY jobid LIMIT 1 OFFSET ${PROCESS};")
if [[ -z "$JOBID" ]]; then
    echo "Error: could not determine JOBID for PROCESS=${PROCESS}. The project database may be empty or corrupt." >&2
    exit 1
fi
sqlite3 -noheader -cmd ".mode list" project.db "SELECT cfg FROM configuration WHERE jobid=${JOBID};" > job_config.toml

# Copy the systematics TOML file
if ! copy_to_local_with_retry "$PROJECT/systematics.toml" "systematics.toml" "systematics configuration"; then
    exit 1
fi

# Copy the input data file(s)
mkdir data

# Extract all paths
#full_paths=$(grep -E '"/(pnfs|exp)' job_config.toml | grep -o '"[^"]*"' | sed 's/"//g')
full_paths=$(grep '"/pnfs' job_config.toml | grep -o '"[^"]*"' | sed 's/"//g')
echo "Found $(echo "$full_paths" | wc -l) input files to copy."

# Copy input files
mkdir -p data
for p in $full_paths; do
    echo "Copying input file: $p"
    ifdh cp "$p" data/
done
ls -lrth data/

# Modify the job_config.toml to use local paths
for p in $full_paths; do
    b=$(basename "$p")
    sed -i "s#\"$p\"#\"data/$b\"#g" job_config.toml
done

# Dump some info for debugging
cat job_config.toml
ls -lrth .
ls -lrth data/

#######################################################################
# Run the analysis
#######################################################################

# Run medulla (selection)
./selection/medulla job_config.toml

MEDULLA_STATUS=$?
if [[ $MEDULLA_STATUS -ne 0 ]]; then
    echo "Error: medulla selection failed with exit status ${MEDULLA_STATUS}; output will not be copied." >&2
    exit "$MEDULLA_STATUS"
fi

if [[ ! -f output.root ]]; then
    echo "Error: medulla selection did not create output.root." >&2
    exit 1
fi

OUTPUT_SIZE=$(stat -c%s output.root)
if [[ $OUTPUT_SIZE -lt 1024 ]]; then
    echo "Error: output.root is only ${OUTPUT_SIZE} bytes; refusing to copy a stub output." >&2
    exit 1
fi


ls -lrth



# Run medulla (systematics)
./systematics/run_systematics systematics.toml
SYSTEMATICS_STATUS=$?
if [[ $SYSTEMATICS_STATUS -ne 0 ]]; then
    echo "Error: systematics failed with exit status ${SYSTEMATICS_STATUS}; outputs will not be copied." >&2
    exit "$SYSTEMATICS_STATUS"
fi

if [[ ! -f output_sys.root ]]; then
    echo "Error: systematics did not create output_sys.root." >&2
    exit 1
fi

SYSTEMATICS_SIZE=$(stat -c%s output_sys.root)
if [[ $SYSTEMATICS_SIZE -lt 1024 ]]; then
    echo "Error: output_sys.root is only ${SYSTEMATICS_SIZE} bytes; refusing to copy a stub output." >&2
    exit 1
fi

ls -lrth

# Copy the systematics output first. The raw output is copied last because
# project status uses output_jobid*.root as the job-completion marker.
printf -v SYSTNAME "output_systematics_jobid%04d.root" "$JOBID"
if ! copy_to_remote_with_retry output_sys.root "$PROJECT/output/$SYSTNAME" "systematics output"; then
    exit 1
fi

printf -v RAWNAME "output_jobid%04d.root" "$JOBID"
if ! copy_to_remote_with_retry output.root "$PROJECT/output/$RAWNAME" "selection output"; then
    exit 1
fi
