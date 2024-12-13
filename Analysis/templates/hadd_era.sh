#!/bin/sh
# Usage:
#    hadd_era.sh

echo "Hadd YEARNAME:"
hadd Analyzed_Data_YEARNAME.root stream_*/Analyzed_Data_*.root
echo "Done YEARNAME!"
