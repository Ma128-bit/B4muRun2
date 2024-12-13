# !/bin/sh
# Usage:
#    prepare_and_submit_ALL.sh <Year> <Delta>

helpstring="Usage:
prepare_and_submit_ALL.sh [Year] [Analysis_type] [Delta]"

year=$1
analysis_type=$2
delta=$3

# Check inputs
if [ -z ${3+x} ]; then
    echo -e ${helpstring}
    return
fi

declare -a Era_2022=("C" "D-v1" "D-v2" "E" "F" "G")
declare -a Era_2023=("B" "C-v1" "C-v2" "C-v3" "C-v4" "D-v1" "D-v2")
declare -a Era_2024=("B" "C" "D" "E-v1" "E-v2" "F" "G" "H" "I-v1" "I-v2")

if [ "${year}" == "2022" ]; then
  eras=("${Era_2022[@]}")
elif [ "${year}" == "2023" ]; then
  eras=("${Era_2023[@]}")
elif [ "${year}" == "2024" ]; then
  eras=("${Era_2024[@]}")
else
  echo "Error: The year is incorrect."
  return
fi

for e in "${eras[@]}"; do
    source prepare_condor.sh $e $year $analysis_type $delta
done

echo "End preparation ... beginning submission"
sleep 1

cd $analysis_type
for e in "${eras[@]}"; do
    cd "${year}_era${e}"
    source submit_era.sh
    echo ""
    cd ..
    sleep 1
done
cd ..
