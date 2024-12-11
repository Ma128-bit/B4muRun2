#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Usage: $0 <type> <year>"
  return
fi

type=$1
year=$2

directory="$PWD"
#echo "pwd: $directory"
home_dir=$(dirname "$(dirname "$directory")/CrabSubmission")
#echo "Home dir: $home_dir"
path_to_skim_file="${home_dir}/SkimTools/test"

declare -a Norm_2017=("Charmonium/Run2017B-UL2017_MiniAODv2-v1" "Charmonium/Run2017C-UL2017_MiniAODv2-v1" "Charmonium/Run2017D-UL2017_MiniAODv2-v1" "Charmonium/Run2017E-UL2017_MiniAODv2-v1" "Charmonium/Run2017F-UL2017_MiniAODv2-v1")
declare -a Era_2017=("B" "C" "D" "E" "F")
declare -a Norm_2018=("Charmonium/Run2018A-UL2018_MiniAODv2-v1" "Charmonium/Run2018B-UL2018_MiniAODv2-v1" "Charmonium/Run2018C-UL2018_MiniAODv2-v1" "Charmonium/Run2018D-UL2018_MiniAODv2-v1")
declare -a Era_2018=("A" "B" "C" "D")

declare -a MC22_B4mu_pre=("/BdTo4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v3/MINIAODSIM" "/Bs0To4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v3/MINIAODSIM" "/BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM")
declare -a MC22_B4mu_post=("/BdTo4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v3/MINIAODSIM" "/Bs0To4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v3/MINIAODSIM" "/BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM")

declare -a MC23_B4mu_pre=("/BdTo4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/MINIAODSIM" "/Bs0To4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/MINIAODSIM" "/BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/MINIAODSIM")
declare -a MC23_B4mu_post=("/BdTo4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3/MINIAODSIM" "/Bs0To4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3/MINIAODSIM" "/BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v2/MINIAODSIM")

declare -a B4mu_MC_label=("Bd" "Bs" "BsJPsiPhi")

if [ "${year}" == "2017" ]; then
    case "$type" in
      Norm)
        Data_ID=("${Norm_2017[@]}")
        eras=("${Era_2017[@]}")
        golden_json="Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
        ;;
      MC)
        globaltag="130X_mcRun3_2022_realistic_v5"
        datasets=("${MC22_B4mu_pre[@]}")
        label=("${B4mu_MC_label[@]}")
        input_type="global"
        ;;
      *)
        echo "Error: The type is incorrect."
        return
        ;;
    esac
elif [ "${year}" == "2018" ]; then
    case "$type" in
      Norm)
        Data_ID=("${Norm_2018[@]}")
        eras=("${Era_2018[@]}")
        golden_json="Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
        ;;
      MC)
        globaltag="130X_mcRun3_2022_realistic_v5"
        datasets=("${MC22_B4mu_pre[@]}")
        label=("${B4mu_MC_label[@]}")
        input_type="global"
        ;;
      *)
        echo "Error: The type is incorrect."
        return
        ;;
    esac
else
    echo "Error: The year is incorrect."
    return
fi

#voms-proxy-init --valid 192:00 --voms cms

if [[ "$type" != *"MC"* ]]; then
    mkdir -p "${year}_type${type}"
    echo "Data ${year} - type ${type} is selected"
    path="${directory}/${year}_type${type}/PatAndTree_cfg.py"
    cp "${path_to_skim_file}/run_Data2022_PatAndTree_cfg.py" "$path"
    #sed -i "s#106X_dataRun2_v37#${globaltag}#g" "${year}_era${era}/PatAndTree_cfg.py"
    cp templates/report.sh "${year}_type${type}/report.sh"
    cp templates/status.sh "${year}_type${type}/status.sh"
    cp templates/resubmit.sh "${year}_type${type}/resubmit.sh"
    cd "${year}_type${type}"
    sed -i "s#YEAR#${year}#g" *.sh
    sed -i "s#ERANAME#${type}#g" *.sh
    cd ..
    cp templates/submit.sh "${year}_type${type}/submit.sh"
    for i in $(seq 0 $((${#eras[@]} - 1))); do
        cp templates/CRAB_template.py "${year}_type${type}/CRAB_stream_${eras[${i}]}.py"
        sed -i "s#YEAR#${year}#g" "${year}_type${type}/CRAB_stream_${eras[${i}]}.py"
        sed -i "s#ERANAME#${type}#g" "${year}_type${type}/CRAB_stream_${eras[${i}]}.py"
        sed -i "s#NUMBER#${i}#g" "${year}_type${type}/CRAB_stream_${eras[${i}]}.py"
        sed -i "s#DATASET_ID#${Data_ID[${i}]}#g" "${year}_type${type}/CRAB_stream_${eras[${i}]}.py"
        sed -i "s#FILE_TO_SUBMIT_PATH#${path}#g" "${year}_type${type}/CRAB_stream_${eras[${i}]}.py"
        sed -i "s#GOLDEN_JSON_PATH#${golden_json}#g" "${year}_type${type}/CRAB_stream_${eras[${i}]}.py"
        cd "${year}_type${type}"
        crab submit -c "CRAB_stream_${eras[${i}]}.py"
        cd ..
        echo "Stream $i submitted!"
        sleep 3
    done
else
    mkdir -p "${year}_${era}"
    echo "${era} - ${year} is selected"
    path="${directory}/${year}_${era}/PatAndTree_cfg.py"
    cp "${path_to_skim_file}/run_MC2022_PatAndTree_cfg.py" "$path"
    sed -i "s#130X_mcRun3_2022_realistic_postEE_v6#${globaltag}#g" "${year}_${era}/PatAndTree_cfg.py"
    j=0
    for i in "${datasets[@]}"; do
        cp templates/CRAB_template_MC.py "${year}_${era}/CRAB_MC_${label[${j}]}.py"
        sed -i "s#YEAR#${year}#g" "${year}_${era}/CRAB_MC_${label[${j}]}.py"
        sed -i "s#ERANAME#${era}#g" "${year}_${era}/CRAB_MC_${label[${j}]}.py"
        sed -i "s#MC_DATASET#${i}#g" "${year}_${era}/CRAB_MC_${label[${j}]}.py"
        sed -i "s#FILE_TO_SUBMIT_PATH#${path}#g" "${year}_${era}/CRAB_MC_${label[${j}]}.py"
        sed -i "s#INPUT_TYPE#${input_type}#g" "${year}_${era}/CRAB_MC_${label[${j}]}.py"
        sed -i "s#B_TYPE#${label[${j}]}#g" "${year}_${era}/CRAB_MC_${label[${j}]}.py"
        cd "${year}_${era}"
        #crab submit -c "CRAB_MC_${label[${j}]}.py"
        cd ..
        echo "${era} - $j submitted!"
        ((j++))
        sleep 3
    done

fi




