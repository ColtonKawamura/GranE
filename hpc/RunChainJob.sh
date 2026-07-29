#!/bin/bash
# Usage: bash RunChainJob.sh RunChainGPU.sh {wc -l jobScript}
# Splits total jobs into chunks of 1000 and chains them
# change the batch lext to the END based on memory andCPU limit

SCRIPT=$1
TOTAL=$2
CHUNK=1000

PREV_JOB=""
START=1

while [ $START -le $TOTAL ]; do
    END=$((START + CHUNK - 1))
    if [ $END -gt $TOTAL ]; then END=$TOTAL; fi

    COUNT=$((END - START + 1))
    OFFSET=$((START - 1))

    if [ -z "$PREV_JOB" ]; then
        JOB=$(sbatch --array=1-${COUNT}%50 --export=ALL,OFFSET=$OFFSET $SCRIPT | awk '{print $NF}')
    else
        JOB=$(sbatch --array=1-${COUNT}%50 --export=ALL,OFFSET=$OFFSET --dependency=afterok:$PREV_JOB $SCRIPT | awk '{print $NF}')
    fi

    echo "Submitted jobs ${START}-${END} as job ID $JOB"
    PREV_JOB=$JOB
    START=$((END + 1))
done

