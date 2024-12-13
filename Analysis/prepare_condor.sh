# !/bin/sh
# Usage:
#    prepare_condor.sh <Year> <Analysis_type> <Delta>

helpstring="Usage:
prepare_condor.sh [Year] [Analysis_type] [Delta]"
year=$1
Analysis_type=$2
delta=$3

# Check inputs
if [ -z ${3+x} ]; then
    echo -e ${helpstring}
    return
fi

file_directory="/lustre/cms/store/user/mbuonsan"

declare -a Data_2017=("SkimB4Mu_2017eraB_Norm_Mini/241211_175616" "SkimB4Mu_2017eraC_Norm_Mini/241211_175626" "SkimB4Mu_2017eraD_Norm_Mini/241211_175637" "SkimB4Mu_2017eraE_Norm_Mini/241211_175649" "SkimB4Mu_2017eraF_Norm_Mini/241211_175658")

declare -a Data_2018=("SkimB4Mu_2018eraA_Norm_Mini/241211_173436" "SkimB4Mu_2018eraB_Norm_Mini/241211_173446" "SkimB4Mu_2018eraC_Norm_Mini/241211_173455" "SkimB4Mu_2018eraD_Norm_Mini/241211_173504")


if [ "${Analysis_type}" == "Norm" ]; then
    datasets="Charmonium"
    if [ "${year}" == "2017" ]; then
        files=("${Data_2017[@]}")
    elif [ "${year}" == "2018" ]; then
        files=("${Data_2018[@]}")
    else
        echo "Error: The year is incorrect."
        return
    fi
elif [ "${Analysis_type}" == "B4mu" ]; then
    datasets="DoubleMuonLowMass"
    if [ "${year}" == "2017" ]; then
        files=("${Data_2017[@]}")
    else
        echo "Error: The year is incorrect."
        return
    fi
fi

home_directory="$PWD"

if [[ "$Analysis_type" != *"MC"* ]]; then
    if [ ! -d "${home_directory}/${Analysis_type}/${year}" ]; then
        mkdir -p "${home_directory}/${Analysis_type}/${year}"
    fi
    echo "Data ${year} is selected"
    cp templates/submit_era.sh "${home_directory}/${Analysis_type}/${year}"
    cp templates/hadd_era.sh "${home_directory}/${Analysis_type}/${year}"
    sed -i "s#YEARNAME#${year}#g" "${home_directory}/${Analysis_type}/${year}/submit_era.sh"
    sed -i "s#YEARNAME#${year}#g" "${home_directory}/${Analysis_type}/${year}/hadd_era.sh"

    for i in $(seq 0 $((${#files[@]} - 1))); do
        if [ ! -d "${home_directory}/${Analysis_type}/${year}/stream_${i}" ]; then
            mkdir -p "${home_directory}/${Analysis_type}/${year}/stream_${i}"
            mkdir -p "${home_directory}/${Analysis_type}/${year}/stream_${i}/log"
        fi

        cp templates/submit.condor "${home_directory}/${Analysis_type}/${year}/stream_${i}"
        ndir=$(ls "${file_directory}/${datasets}/${files[${i}]}/" | wc -l)
        tot=$(find "${file_directory}/${datasets}/${files[${i}]}" -type f | wc -l)
        number_of_splits=$(((${tot} / ${delta}) + 1))
        echo "queue ${number_of_splits}" >> "${home_directory}/${Analysis_type}/${year}/stream_${i}/submit.condor"
        sed -i "s#PATH#${home_directory}/${Analysis_type}/${year}/stream_${i}#g" "${home_directory}/${Analysis_type}/${year}/stream_${i}/submit.condor"
        chmod a+x "${home_directory}/${Analysis_type}/${year}/stream_${i}/submit.condor"
        
        cp templates/launch_analysis.sh "${home_directory}/${Analysis_type}/${year}/stream_${i}"
        sed -i "s#PATH#${home_directory}#g" "${home_directory}/${Analysis_type}/${year}/stream_${i}/launch_analysis.sh"
        sed -i "s#DELTAVAL#${delta}#g" "${home_directory}/${Analysis_type}/${year}/stream_${i}/launch_analysis.sh"
        sed -i "s#INPUT_DIR#${file_directory}/${datasets}/${files[${i}]}/#g" "${home_directory}/${Analysis_type}/${year}/stream_${i}/launch_analysis.sh"
        sed -i "s#OUTPUT_DIR#${home_directory}/${Analysis_type}/${year}/stream_${i}#g" "${home_directory}/${Analysis_type}/${year}/stream_${i}/launch_analysis.sh"
        sed -i "s#TRUEFALSE#0#g" "${home_directory}/${Analysis_type}/${year}/stream_${i}/launch_analysis.sh"
        sed -i "s#ANALYSISTYPE#${Analysis_type}#g" "${home_directory}/${Analysis_type}/${year}/stream_${i}/launch_analysis.sh"
        chmod a+x "${home_directory}/${Analysis_type}/${year}/stream_${i}/launch_analysis.sh"
        
        echo -n "."
        sleep 1
    done
else
    if [ ! -d "${home_directory}/${Analysis_type}/${year}_${era}" ]; then
        mkdir -p "${home_directory}/${Analysis_type}/${year}_${era}"
    fi
    echo "Data ${year} - ${era} is selected"
    j=0
    for i in "${datasets[@]}"; do
        if [ ! -d "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}" ]; then
            mkdir -p "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}"
            mkdir -p "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/log"
        fi
        cp templates/submit.condor "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}"
        ndir=$(ls "${file_directory}/${i}/" | wc -l)
        tot=0
        for k in $(seq 0 $((ndir - 1))); do
            nfiles=$(ls "${file_directory}/${i}/000${k}/" | wc -l)
            tot=$((tot + nfiles))
        done
        number_of_splits=$(((${tot} / ${delta}) + 1))
        echo "queue ${number_of_splits}" >> "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/submit.condor"
        sed -i "s#PATH#${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}#g" "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/submit.condor"
        chmod a+x "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/submit.condor"
        
        cp templates/launch_analysis.sh "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}"
        sed -i "s#PATH#${home_directory}#g" "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/launch_analysis.sh"
        sed -i "s#DELTAVAL#${delta}#g" "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/launch_analysis.sh"
        sed -i "s#INPUT_DIR#${file_directory}/${i}#g" "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/launch_analysis.sh"
        sed -i "s#OUTPUT_DIR#${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}#g" "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/launch_analysis.sh"
        if [[ "$i" == *"BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen"* ]]; then
            sed -i "s#TRUEFALSE#2#g" "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/launch_analysis.sh"
        else
            sed -i "s#TRUEFALSE#1#g" "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/launch_analysis.sh"
        fi
        sed -i "s#ANALYSISTYPE#${Analysis_type}#g" "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/launch_analysis.sh"
        chmod a+x "${home_directory}/${Analysis_type}/${year}_${era}/${label[${j}]}/launch_analysis.sh"
        
        echo -n "."
        ((j++))
        sleep 1
    done
fi
echo " Done!"


