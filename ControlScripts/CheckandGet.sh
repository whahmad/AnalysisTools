#! /bin/bash
# This script retrieves jobs from the grid and returns the number of jobs check and recieved
if [ "${1}"  == "--help" ] || [ "${1}" == "" ]; then
    echo "Options for running this script:"
    echo "--help                    Prints this message"
    echo "--get                     Retrieve completed jobs"
    echo "--status                  Show status of submitted jobs"
fi
if [ "${1}"  == "--status" ]; then
    cat jobs_submitted | awk '{system("glite-ce-job-status " $1 )}' > status
    echo '# of jobs with status DONE-OK'
    cat status | grep -c DONE-OK
    echo '# of jobs with status DONE-FAILED'
    cat status | grep -c DONE-FAILED
    echo '# of jobs with status REALLY-RUNNING'
    cat status | grep -c REALLY-RUNNING
    echo '# of jobs with status RUNNING'
    cat status | grep -c RUNNING
    echo '# of jobs with status IDLE'
    cat status | grep -c IDLE
    echo '# of jobs with status CANCELLED'
    cat status | grep -c CANCELLED
    echo '# of jobs with status REGISTERED'
    cat status | grep -c REGISTERED
    rm status
fi
if [ "${1}"  == "--get" ]; then
    touch jobs_complete
    touch jobs_get

    cat jobs_submitted | awk '{system("glite-ce-job-status " $1 " | grep -i done | wc -l |awk \x27{sum++}END{if($1 == 1) print  \" glite-ce-job-output --dir "$2 " " $1 " ; mv " $2 "/*_CREAM*/* "$2 "; rmdir " $2 "/*_CREAM* ; echo " $1 " " $2 "  >> jobs_complete; source " $2 "/GRIDRetrieve.sh  \"}\x27 >> jobs_get ")}'

    source jobs_get 
    rm jobs_get
    cat jobs_complete | awk '{system("grep -v \"" $1 " " $2 "\" jobs_submitted > junk ; cp junk jobs_submitted; rm junk;")}'
    echo 'Number of jobs produced but not complete and downloaded:'  
    cat  jobs_submitted  | wc -l
    echo 'Number of jobs Complete and Downloaded:'
    cat  jobs_complete  | wc -l
fi

