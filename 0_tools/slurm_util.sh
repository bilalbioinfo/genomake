wait_for_slurm_job() {
    local job_id="$1"

    printf "\t\tWaiting for job $job_id to finish...\n"

    while true; do
        local status=$(sacct -j "$job_id" --format=State --noheader | awk '{print $1}' | head -n 1)

        if [[ "$status" =~ COMPLETED|FAILED|CANCELLED|TIMEOUT ]]; then
            break
        fi

        sleep 10
    done

    if [[ "$status" == "COMPLETED" ]]; then
        printf "\t\tJob $job_id finished with status: $status\n"
    else
        printf "\t\tERROR: Job $job_id failed (status: $status)\n"
        exit 1
    fi
}