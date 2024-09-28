HISTOUTPUT=$PWD/output

if [ $# -eq 0 ]; then
    hadd ${HISTOUTPUT}/output_all.root ${HISTOUTPUT}/output_*.root 
else
    rm -f myhadd.list
    # move all existing files to old file

    directory_suffix=1
    base_directory=${HISTOUTPUT}/unmerged_files

    # find the newest directory name that exist
    while true; do
        target_directory="${base_directory}${directory_suffix}"

        if [ ! -d "${target_directory}" ]; then
            mkdir -p "${target_directory}"
            echo "All unmerged root files are moved to " ${target_directory}
            mv ${HISTOUTPUT}/output_*.root "${target_directory}"
            break
        else
            ((directory_suffix++))
        fi
    done

    ls "${target_directory}"/output_*.root > myhadd.list
    cp myhadd.C ${target_directory}

    root -b -l -q 'myhadd.cpp("'${target_directory}/'",0,"'${HISTOUTPUT}'",'$1')'
fi
