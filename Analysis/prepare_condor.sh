# !/bin/sh
# Usage:
#    prepare_condor.sh <Era> <Year> <Analysis_type> <Delta>

helpstring="Usage:
prepare_condor.sh [Era] [Year] [Analysis_type] [Delta]"
era=$1
year=$2
Analysis_type=$3
delta=$4

# Check inputs
if [ -z ${4+x} ]; then
    echo -e ${helpstring}
    return
fi

if [ "${Analysis_type}" == "B2mu2K" ]; then
    Ana_temp="B2Mu2K"
elif [ "${Analysis_type}" == "B2muKpi" ]; then
    Ana_temp="B2Mu2K"
elif [ "${Analysis_type}" == "B4mu" ]; then
    Ana_temp="B4Mu"
else
    echo "Error: The Analysis_type is incorrect."
    return
fi

file_directory="/lustre/cms/store/user/mbuonsan"

declare -a C_2022=("240731_160611" "240731_160630" "240731_160644" "240731_162051" "240731_160721" "240731_160745" "240731_160804" "240731_160818")
declare -a D_v1_2022=("240731_160834" "240731_160859" "240731_160917" "240731_160933" "240731_160948" "240731_161012" "240731_161027" "240731_161042")
declare -a D_v2_2022=("240731_161059" "240731_161122" "240731_161136" "240731_161149" "240731_161205" "240731_161223" "240731_161240" "240731_161253")
declare -a E_2022=("240731_161312" "240731_161332" "240731_161352" "240731_161411" "240731_161426" "240731_161440" "240731_161506" "240731_161522")
declare -a F_2022=("240731_153720" "240731_153734" "240731_153753" "240731_153810" "240731_153823" "240731_153837" "240731_153854" "240731_153916")
declare -a G_2022=("240731_153931" "240731_153945" "240731_153959" "240731_154014" "240731_154037" "240731_154052" "240731_154111" "240731_154127")

declare -a B_2023=("240731_154536" "240731_154557" "240731_154619" "240731_154643" "240731_154707" "240731_154721" "240731_154741" "240731_154756")
declare -a C_v1_2023=("240731_154810" "240731_154827" "240731_154847" "240731_154902" "240731_154915" "240731_154931" "240731_154944" "240731_155009")
declare -a C_v2_2023=("240731_155026" "240731_155044" "240731_155100" "240731_155117" "240731_155133" "240731_155149" "240731_155206" "240731_155229")
declare -a C_v3_2023=("240731_155245" "240731_155304" "240731_155320" "240731_155344" "240731_155401" "240731_155427" "240731_155454" "240731_155509")
declare -a C_v4_2023=("240731_155534" "240731_155557" "240731_155613" "240731_155630" "240731_155654" "240731_155714" "240731_155729" "240731_155745")
declare -a D_v1_2023=("240731_155815" "240731_155830" "240731_155844" "240731_155858" "240731_155923" "240731_155938" "240731_155952" "240731_160008")
declare -a D_v2_2023=("240731_160031" "240731_160049" "240731_160104" "240731_160117" "240731_160138" "240731_160158" "240731_160213" "240731_160227")

declare -a B_2024=("241122_174414" "241125_092457" "241125_092508" "241125_092520" "241125_092532" "241125_092546" "241125_092558" "241125_092610")
declare -a C_2024=("241125_092710" "241125_092722" "241125_092732" "241125_092742" "241125_092757" "241125_092807" "241125_092818" "241125_092827")
declare -a D_2024=("241125_092838" "241125_092847" "241125_092858" "241125_092910" "241125_092920" "241125_092931" "241125_092942" "241125_092951")
declare -a E_v1_2024=("241125_093003" "241125_093013" "241125_093023" "241125_093033" "241125_093044" "241125_093054" "241125_093103" "241125_093113")
declare -a E_v2_2024=("241125_093126" "241125_093136" "241125_093145" "241125_093155" "241125_093206" "241125_093217" "241125_093228" "241125_093239")
declare -a F_2024=("241125_093252" "241125_093302" "241125_093312" "241125_093322" "241125_093332" "241125_093342" "241125_093351" "241125_093401")
declare -a G_2024=("241125_093412" "241125_093422" "241125_093432" "241125_093442" "241125_093452" "241125_093502" "241125_093513" "241125_093522")
declare -a H_2024=("241125_093534" "241125_093544" "241125_093553" "241125_093606" "241125_093616" "241125_093626" "241125_093635" "241125_093644")
declare -a I_v1_2024=("241125_093658" "241125_093711" "241125_093723" "241125_093737" "241125_093747" "241125_093757" "241125_093807" "241125_093817")
declare -a I_v2_2024=("241125_093828" "241125_093839" "241125_093848" "241125_093858" "241125_093908" "241125_093918" "241125_093927" "241125_093937")

declare -a MC22_B4mu_pre=("BdTo4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimB4Mu_2022_MC_pre_Bd_Mini/240731_153207" "Bs0To4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimB4Mu_2022_MC_pre_Bs_Mini/240731_153225" "BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/SkimB4Mu_2022_MC_pre_BsJPsiPhi_Mini/240731_153250")
declare -a MC22_B4mu_post=("BdTo4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimB4Mu_2022_MC_post_Bd_Mini/240731_153311" "Bs0To4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimB4Mu_2022_MC_post_Bs_Mini/240731_153325" "BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/SkimB4Mu_2022_MC_post_BsJPsiPhi_Mini/240731_153342")

declare -a MC23_B4mu_pre=("BdTo4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimB4Mu_2023_MC_pre_Bd_Mini/240731_153505" "Bs0To4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimB4Mu_2023_MC_pre_Bs_Mini/240731_153523" "BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/SkimB4Mu_2023_MC_pre_BsJPsiPhi_Mini/240731_153546")
declare -a MC23_B4mu_post=("BdTo4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimB4Mu_2023_MC_post_Bd_Mini/240731_153610" "Bs0To4Mu_FourMuonFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimB4Mu_2023_MC_post_Bs_Mini/240731_153630" "BsToJpsiPhi_JMM_PhiMM_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/SkimB4Mu_2023_MC_post_BsJPsiPhi_Mini/240731_153650")

declare -a B4mu_MC_label=("Bd_4mu" "Bs_4mu" "BsJPsiPhi")

declare -a Control_C_2022=("241122_165954" "241122_171707" "241122_173419" "241124_100840" "241124_102559" "241124_103315" "241124_103530" "241124_103535")
declare -a Control_D_v1_2022=("241124_104607" "241124_110319" "241124_111259" "241124_111248" "241124_111250" "241124_105402" "241124_105414" "241124_105433")
declare -a Control_D_v2_2022=("241124_113223" "241124_113229" "241124_113235" "241124_113235" "241124_113247" "241124_113257" "241124_113259" "241124_113309")
declare -a Control_E_2022=("241125_194542" "241124_115003" "241124_115009" "241124_115015" "241124_115045" "241124_115045" "241124_115049" "241124_115056")
declare -a Control_F_2022=("241122_122123" "241122_122131" "241122_122138" "241122_122145" "241122_122152" "241122_122159" "241122_122207" "241122_122214")
declare -a Control_G_2022=("241122_122228" "241122_122235" "241122_122242" "241122_122249" "241122_122256" "241122_122304" "241122_122311" "241122_122317")

declare -a Control_B_2023=("241122_122624" "241122_122631" "241122_122638" "241122_122644" "241122_122651" "241122_122657" "241122_122704" "241122_122711")
declare -a Control_C_v1_2023=("241122_122721" "241122_122728" "241122_122735" "241122_122741" "241122_122748" "241122_122755" "241122_122802" "241122_122808")
declare -a Control_C_v2_2023=("241122_122817" "241122_122823" "241122_122830" "241122_122837" "241122_122843" "241122_122850" "241122_122857" "241122_122904")
declare -a Control_C_v3_2023=("241122_122913" "241122_122920" "241122_122927" "241122_122934" "241122_122940" "241122_122947" "241122_122954" "241122_123001")
declare -a Control_C_v4_2023=("241122_123009" "241122_123016" "241122_123023" "241122_123029" "241122_123036" "241122_123043" "241122_123050" "241122_123057")
declare -a Control_D_v1_2023=("241122_123107" "241122_123114" "241122_123121" "241122_123128" "241122_123135" "241122_123142" "241122_123149" "241122_123156")
declare -a Control_D_v2_2023=("241122_123205" "241122_123212" "241122_123220" "241122_123227" "241122_123234" "241122_123241" "241122_123256" "241122_123304")

declare -a MC22_B2mu2trk_pre=("BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2022_MC_pre_Bs2mu2K_Mini/241122_122443" "BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2022_MC_pre_Bd2muKpi_Mini/241122_122449")
declare -a MC22_B2mu2trk_post=("BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2022_MC_post_Bs2mu2K_Mini/241122_122459" "BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2022_MC_post_Bd2muKpi_Mini/241122_122505")

declare -a MC23_B2mu2trk_pre=("BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2023_MC_pre_Bs2mu2K_Mini/241122_122536" "BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2023_MC_pre_Bd2muKpi_Mini/241122_122543")
declare -a MC23_B2mu2trk_post=("BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2023_MC_post_Bs2mu2K_Mini/241122_122554" "BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2023_MC_post_Bd2muKpi_Mini/241122_122600")

declare -a B2mu2trk_MC_label=("B2mu2K" "B2muKpi")

declare -a MC22_B2mu2K_pre_new=("BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2022_MC_pre_new_Bs2mu2K_Mini/241203_142000" "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2022_MC_pre_new_ext_Bs2mu2K_Mini/241203_142039")
declare -a MC22_B2mu2K_post_new=("BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2022_MC_post_new_Bs2mu2K_Mini/241203_142018" "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2022_MC_post_new_ext_Bs2mu2K_Mini/241203_142052")
declare -a B2mu2K22_MC_label=("B2mu2K_new" "B2mu2K_new_ext")
declare -a MC23_B2mu2K_pre_new=("BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2023_MC_pre_new_Bs2mu2K_Mini/241203_142109")
declare -a MC23_B2mu2K_post_new=("BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/SkimB2Mu2K_2023_MC_post_new_Bs2mu2K_Mini/241203_142129")
declare -a B2mu2K23_MC_label=("B2mu2K_new")

if [ "${Ana_temp}" == "B4Mu" ]; then
    if [ "${year}" == "2022" ]; then
        case "$era" in
          C)
            datasets=("${C_2022[@]}")
            ;;
          D-v1)
            datasets=("${D_v1_2022[@]}")
            ;;
          D-v2)
            datasets=("${D_v2_2022[@]}")
            ;;
          E)
            datasets=("${E_2022[@]}")
            ;;
          F)
            datasets=("${F_2022[@]}")
            ;;
          G)
            datasets=("${G_2022[@]}")
            ;;
          MC_pre)
            datasets=("${MC22_B4mu_pre[@]}")
            label=("${B4mu_MC_label[@]}")
            ;;
          MC_post)
            datasets=("${MC22_B4mu_post[@]}")
            label=("${B4mu_MC_label[@]}")
            ;;
          *)
            echo "Error: The era is incorrect."
            return
            ;;
        esac
    elif [ "${year}" == "2023" ]; then
        case "$era" in
          B)
            datasets=("${B_2023[@]}")
            ;;
          C-v1)
            datasets=("${C_v1_2023[@]}")
            ;;
          C-v2)
            datasets=("${C_v2_2023[@]}")
            ;;
          C-v3)
            datasets=("${C_v3_2023[@]}")
            ;;
          C-v4)
            datasets=("${C_v4_2023[@]}")
            ;;
          D-v1)
            datasets=("${D_v1_2023[@]}")
            ;;
          D-v2)
            datasets=("${D_v2_2023[@]}")
            ;;
          MC_pre)
            datasets=("${MC23_B4mu_pre[@]}")
            label=("${B4mu_MC_label[@]}")
            ;;
          MC_post)
            datasets=("${MC23_B4mu_post[@]}")
            label=("${B4mu_MC_label[@]}")
            ;;
          *)
            echo "Error: The era is incorrect."
            return
            ;;
        esac
    elif [ "${year}" == "2024" ]; then
        case "$era" in
          B)
            datasets=("${B_2024[@]}")
            ;;
          C)
            datasets=("${C_2024[@]}")
            ;;
          D)
            datasets=("${D_2024[@]}")
            ;;
          E-v1)
            datasets=("${E_v1_2024[@]}")
            ;;
          E-v2)
            datasets=("${E_v2_2024[@]}")
            ;;
          F)
            datasets=("${F_2024[@]}")
            ;;
          G)
            datasets=("${G_2024[@]}")
            ;;
          H)
            datasets=("${H_2024[@]}")
            ;;
          I-v1)
            datasets=("${I_v1_2024[@]}")
            ;;
          I-v2)
            datasets=("${I_v2_2024[@]}")
            ;;
          *)
            echo "Error: The era is incorrect."
            return
            ;;
        esac
    else
        echo "Error: The year is incorrect."
        return
    fi
elif [ "${Ana_temp}" == "B2Mu2K" ]; then
    if [ "${year}" == "2022" ]; then
        case "$era" in
          C)
            datasets=("${Control_C_2022[@]}")
            ;;
          D-v1)
            datasets=("${Control_D_v1_2022[@]}")
            ;;
          D-v2)
            datasets=("${Control_D_v2_2022[@]}")
            ;;
          E)
            datasets=("${Control_E_2022[@]}")
            ;;
          F)
            datasets=("${Control_F_2022[@]}")
            ;;
          G)
            datasets=("${Control_G_2022[@]}")
            ;;
          MC_pre)
            datasets=("${MC22_B2mu2trk_pre[@]}")
            label=("${B2mu2trk_MC_label[@]}")
            ;;
          MC_post)
            datasets=("${MC22_B2mu2trk_post[@]}")
            label=("${B2mu2trk_MC_label[@]}")
            ;;
          MC_pre_new)
            datasets=("${MC22_B2mu2K_pre_new[@]}")
            label=("${B2mu2K22_MC_label[@]}")
            ;;
          MC_post_new)
            datasets=("${MC22_B2mu2K_post_new[@]}")
            label=("${B2mu2K22_MC_label[@]}")
            ;;
          *)
            echo "Error: The era is incorrect."
            return
            ;;
        esac
    elif [ "${year}" == "2023" ]; then
        case "$era" in
          B)
            datasets=("${Control_B_2023[@]}")
            ;;
          C-v1)
            datasets=("${Control_C_v1_2023[@]}")
            ;;
          C-v2)
            datasets=("${Control_C_v2_2023[@]}")
            ;;
          C-v3)
            datasets=("${Control_C_v3_2023[@]}")
            ;;
          C-v4)
            datasets=("${Control_C_v4_2023[@]}")
            ;;
          D-v1)
            datasets=("${Control_D_v1_2023[@]}")
            ;;
          D-v2)
            datasets=("${Control_D_v2_2023[@]}")
            ;;
          MC_pre)
            datasets=("${MC23_B2mu2trk_pre[@]}")
            label=("${B2mu2trk_MC_label[@]}")
            ;;
          MC_post)
            datasets=("${MC23_B2mu2trk_post[@]}")
            label=("${B2mu2trk_MC_label[@]}")
            ;;
          MC_pre_new)
            datasets=("${MC23_B2mu2K_pre_new[@]}")
            label=("${B2mu2K23_MC_label[@]}")
            ;;
          MC_post_new)
            datasets=("${MC23_B2mu2K_post_new[@]}")
            label=("${B2mu2K23_MC_label[@]}")
            ;;
          *)
            echo "Error: The era is incorrect."
            return
            ;;
        esac
    else
        echo "Error: The year is incorrect."
        return
    fi
fi


home_directory="$PWD"

if [[ "$era" != *"MC"* ]]; then
    if [ ! -d "${home_directory}/${Analysis_type}/${year}_era${era}" ]; then
        mkdir -p "${home_directory}/${Analysis_type}/${year}_era${era}"
    fi
    echo "Data ${year} - era ${era} is selected"
    cp templates/submit_era.sh "${home_directory}/${Analysis_type}/${year}_era${era}"
    cp templates/hadd_era.sh "${home_directory}/${Analysis_type}/${year}_era${era}"
    sed -i "s#YEARNAME#${year}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/submit_era.sh"
    sed -i "s#ERANAME#${era}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/submit_era.sh"
    sed -i "s#YEARNAME#${year}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/hadd_era.sh"
    sed -i "s#ERANAME#${era}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/hadd_era.sh"

    for i in {0..7}; do
        if [ ! -d "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}" ]; then
            mkdir -p "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}"
            mkdir -p "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/log"
        fi

        cp templates/submit.condor "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}"
        ndir=$(ls "${file_directory}/ParkingDoubleMuonLowMass${i}/Skim${Ana_temp}_${year}era${era}_stream${i}_Mini/${datasets[${i}]}/" | wc -l)
        tot=$(find "${file_directory}/ParkingDoubleMuonLowMass${i}/Skim${Ana_temp}_${year}era${era}_stream${i}_Mini/${datasets[${i}]}" -type f | wc -l)
        #wtot=$(du -sBG "${file_directory}/ParkingDoubleMuonLowMass${i}/Skim${Ana_temp}_${year}era${era}_stream${i}_Mini/${datasets[${i}]}" | awk '{print $1}' | sed 's/G//')
        #w0=$(echo "$delta * 0.05 / 1" | bc)
        #number_of_splits=$(((${wtot} / ${w0}) + 1))
        #delta=$(((${tot} / ${number_of_splits}) + 1))
        #tot=0
        #for j in $(seq 0 $((ndir - 1))); do
        #    nfiles=$(ls "${file_directory}/ParkingDoubleMuonLowMass${i}/Skim${Ana_temp}_${year}era${era}_stream${i}_Mini/${datasets[${i}]}/000${j}/" | wc -l)
        #    tot=$((tot + nfiles))
        #done
        #echo "nfiles=${tot}"
        number_of_splits=$(((${tot} / ${delta}) + 1))
        echo "queue ${number_of_splits}" >> "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/submit.condor"
        sed -i "s#PATH#${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/submit.condor"
        chmod a+x "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/submit.condor"
        
        cp templates/launch_analysis.sh "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}"
        sed -i "s#PATH#${home_directory}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/launch_analysis.sh"
        sed -i "s#DELTAVAL#${delta}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/launch_analysis.sh"
        sed -i "s#INPUT_DIR#${file_directory}/ParkingDoubleMuonLowMass${i}/Skim${Ana_temp}_${year}era${era}_stream${i}_Mini/${datasets[${i}]}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/launch_analysis.sh"
        sed -i "s#OUTPUT_DIR#${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/launch_analysis.sh"
        sed -i "s#TRUEFALSE#0#g" "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/launch_analysis.sh"
        sed -i "s#ANALYSISTYPE#${Analysis_type}#g" "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/launch_analysis.sh"
        chmod a+x "${home_directory}/${Analysis_type}/${year}_era${era}/stream_${i}/launch_analysis.sh"
        
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


