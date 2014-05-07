#! /bin/bash
# This script gives a job status summary
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

