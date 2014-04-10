#! /bin/bash
cat jobs_submitted | awk '{system("glite-ce-job-purge " $1 " ;")}'

