#!/bin/sh
# Usage:
#    submit_era.sh

echo "Submit YEARNAME:"
for i in {0..MAXNJOB}; do
  echo -n "Stream ${i} "
  condor_submit -name ettore "stream_${i}/submit.condor"
done
echo "Done YEARNAME!"
