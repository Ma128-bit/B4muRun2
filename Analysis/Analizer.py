import os, time
start = time.time()
import argparse
from ROOT import RDataFrame, gROOT, EnableImplicitMT, gInterpreter

gROOT.SetBatch(True)
EnableImplicitMT(16)

gInterpreter.Declare("""
    #include "Utilities.h"
""")

from ROOT import flat2D, flat1D_int, flat1D_double, flat0D_int, flat0D_double, add_int, add_double#, TwoObjMassFit, FourObjMassFit

def list_of_root_files(directory):
    file_root = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".root"):
                rel_path = os.path.relpath(root, directory)
                file_root.append(os.path.join(rel_path, file))
    return file_root

def select_root_files(file_root, i , delta):
    selected_files = []
    for file in file_root:
        temp = ((file.split('_')[1]).split('.'))[0]
        file_ID = int(temp)
        if(file_ID>= i*delta and file_ID < (i+1)*delta):
            selected_files.append(file)
    selected_files = sorted(selected_files)
    #print(selected_files)
    if directory.endswith("/"):
        selected_files = [directory + s for s in selected_files]
    else:
        selected_files = [directory + "/" + s for s in selected_files]
    return(selected_files)

def MuonIDs(rdf, branches, n_muons=4):
    rdf = rdf.Define("Stats","get_stat(mu_index, Muon_isGlobal, Muon_isPF, Muon_isLoose, Muon_isMedium, Muon_isTight, Muon_isSoft, Muon_isTrackerMuon)")
    branches = branches + ["isGlobal", "isPF", "isLoose", "isMedium","isTight", "isSoft", "isTracker"]
    rdf = rdf.Define("isGlobal", flat1D_int(0), ["Stats"])
    rdf = rdf.Define("isPF", flat1D_int(1), ["Stats"])
    rdf = rdf.Define("isLoose", flat1D_int(2), ["Stats"])
    rdf = rdf.Define("isMedium", flat1D_int(3), ["Stats"])
    rdf = rdf.Define("isTight", flat1D_int(4), ["Stats"])
    rdf = rdf.Define("isSoft", flat1D_int(5), ["Stats"])
    rdf = rdf.Define("isTracker", flat1D_int(6), ["Stats"])

    return rdf, branches

def Flat_MuVar(rdf, branches):
    for i in range(1,5):
        ind=str(i)
        for s in ["Pt", "Eta", "Phi"]:
            branches.append("Mu"+ind+"_"+s)
            branches.append("RefTrack"+ind+"_"+s)
            rdf = rdf.Redefine("Mu"+ind+"_"+s,"flattering(Mu"+ind+"_"+s+", Quadruplet_index, \"FLAT MU VAR\")") 
            rdf = rdf.Redefine("RefTrack"+ind+"_"+s,"flattering(RefTrack"+ind+"_"+s+", Quadruplet_index, \"FLAT MU VAR\")") 
    return rdf

def QuadMuVar(rdf, branches, analysis_type):
    # Fix RefittedSV_Mass_err
    quadruplet_related_var = ["Quadruplet_Mass", "FlightDistBS_SV_Significance", "QuadrupletVtx_Chi2", "QuadrupletVtx_NDOF","Quadruplet_Charge", "QuadrupletVtx_x", "QuadrupletVtx_y", "QuadrupletVtx_z", 
                              "RefittedPV_x", "RefittedPV_y", "RefittedPV_z", "Quadruplet_Pt", "Quadruplet_Eta", "Quadruplet_Phi", "FlightDistPVSV", "mu1_pfreliso03", "mu2_pfreliso03", "mu3_pfreliso03", "mu4_pfreliso03", 
                              "mu1_bs_dxy_sig", "mu2_bs_dxy_sig", "mu3_bs_dxy_sig", "mu4_bs_dxy_sig", "vtx_prob", "vtx_ref_prob", "RefittedSV_Chi2", "RefittedSV_nDOF", "RefittedSV_Mass", "RefittedSV_Mass_err", "BS_x", "BS_y", "BS_z"] #FlightDistBS_SV_Significance = lxy_sig        
    vertex_chi2=""
    for i in range(1, 4):
        for j in range(i+1,5):
            vertex_chi2 = vertex_chi2 + ", Vtx"+str(i)+str(j)+"_Chi2"
            quadruplet_related_var.append("Vtx"+str(i)+str(j)+"_Chi2")
            quadruplet_related_var.append("Vtx"+str(i)+str(j)+"_nDOF")
            quadruplet_related_var.append("Vtx"+str(i)+str(j)+"_mass")
            quadruplet_related_var.append("Vtx"+str(i)+str(j)+"_mass_err")
    
    for var in quadruplet_related_var:
        branches.append(var)
        rdf = rdf.Redefine(var,"flattering("+var+", Quadruplet_index, \"FLAT "+var.replace("Vtx","")+" in loop\")") 

    branches.append("Quadruplet_Mass_no_refit") #Not refitted 4mu mass
    rdf = rdf.Define("Quadruplet_Mass_no_refit", "NoRefitMass"+analysis_type+"(MuonPt, Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt, Mu3_Eta, Mu4_Eta, Mu3_Phi, Mu4_Phi, MuonEta, MuonPhi, MuonEnergy)")
    return rdf, vertex_chi2

def MVA_inputs(rdf, branches):
    #cos(θ) angle between B flight direction and 4-muon momentum
    branches.append("Cos3d_PV_SV") #cos3d (not right, use BS)
    rdf = rdf.Define("Cos3d_PV_SV","Cos3D_(QuadrupletVtx_x, QuadrupletVtx_y, QuadrupletVtx_z, RefittedPV_x, RefittedPV_y, RefittedPV_z, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi)")
    branches.append("Cos3d_BS_SV") #cos3d
    rdf = rdf.Define("Cos3d_BS_SV","Cos3D_(QuadrupletVtx_x, QuadrupletVtx_y, QuadrupletVtx_z, BS_x, BS_y, BS_z, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi)")
    
    branches.append("Cos2d_PV_SV") 
    rdf = rdf.Define("Cos2d_PV_SV","Cos2D_(QuadrupletVtx_x, QuadrupletVtx_y, RefittedPV_x, RefittedPV_y, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi)")
    branches.append("Cos2d_BS_SV") #cos2d
    rdf = rdf.Define("Cos2d_BS_SV","Cos2D_(QuadrupletVtx_x, QuadrupletVtx_y, BS_x, BS_y, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi)")

    #∆R max (maximum R distance between any of the 4 muons and the direction of the sum of the 4 muons momenta)
    branches.append("dR_max") #dr
    rdf = rdf.Define("dR_max", "dR_Max(Quadruplet_Eta, Quadruplet_Phi, Mu1_Eta, Mu1_Phi, Mu2_Eta, Mu2_Phi, Mu3_Eta, Mu3_Phi, Mu4_Eta, Mu4_Phi)")
    return rdf

def HLT_quantities(rdf, branches):
    branches = branches + ["HLT_pt_mu1","HLT_pt_mu2","HLT_eta_mu1","HLT_eta_mu2","HLT_phi_mu1","HLT_phi_mu2", "HLT_opposite_charge"]
    rdf = rdf.Define("HLT_pt_mu1", flat0D_double(0), ["MuonPt_HLT"])
    rdf = rdf.Define("HLT_pt_mu2", flat0D_double(1), ["MuonPt_HLT"])
    rdf = rdf.Define("HLT_eta_mu1", flat0D_double(0), ["MuonEta_HLT"])
    rdf = rdf.Define("HLT_eta_mu2", flat0D_double(1), ["MuonEta_HLT"])
    rdf = rdf.Define("HLT_phi_mu1", flat0D_double(0), ["MuonPhi_HLT"])
    rdf = rdf.Define("HLT_phi_mu2", flat0D_double(1), ["MuonPhi_HLT"])
    rdf = rdf.Define("HLT_opposite_charge", "HLT_opposite_charge(MuonPt_HLT, MuonEta_HLT, MuonPhi_HLT, MuonPt, MuonEta, MuonPhi, MuonCharge)")
    return rdf, branches

def DiMuVar(rdf, branches, vertex_chi2):
    #Dimuon masses
    rdf = rdf.Define("Dimuon_index","Dimuon(mu_index, Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt, MuonPt, MuonEta, MuonPhi, MuonCharge)")
    rdf = rdf.Define("Dimuon_mass","DimuonMass(Dimuon_index, MuonPt, MuonEta, MuonPhi, MuonEnergy)")
    #Dimuon dR
    rdf = rdf.Define("Dimuon_dR","DimuondR(Dimuon_index, MuonEta, MuonPhi)")  
    #Dimuon vertex chi2:
    rdf = rdf.Define("Dimuon_chi2","DimuonChi2(mu_index, Dimuon_index, Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt, MuonPt"+ vertex_chi2+")")
    #Flat mass and chi2
    for i in range(2):
        for j in range(2):
            name_mass = "Dimu_OS"+str(i+1)+"_"+str(j+1)
            name_chi2 = "Dimu_OS"+str(i+1)+"_"+str(j+1)+"_chi2"
            branches.append(name_mass)
            branches.append(name_chi2)
            rdf = rdf.Define(name_mass, flat2D(i, j), ["Dimuon_mass"])
            rdf = rdf.Define(name_chi2, flat2D(i, j), ["Dimuon_chi2"])

    branches_toadd = ["Dimu_OS1_dR", "Dimu_OS2_dR", "Quadruplet_Mass_eq", "Dimu_OS_max", "Dimu_OS_min", "isJPsiPhi"]
    for b in branches_toadd:
        branches.append(b)
    
    # Flat Dimuon_dR
    rdf = rdf.Define("Dimu_OS1_dR", flat0D_double(0), ["Dimuon_dR"])
    rdf = rdf.Define("Dimu_OS2_dR", flat0D_double(1), ["Dimuon_dR"])

    #Di muon final Mass
    rdf = rdf.Define("DimuonMassfinal","DimuonMassfinal(Dimu_OS1_1, Dimu_OS1_2, Dimu_OS2_1, Dimu_OS2_2)")
    rdf = rdf.Define("Dimu_OS_max", flat0D_double(0), ["DimuonMassfinal"])
    rdf = rdf.Define("Dimu_OS_min", flat0D_double(1), ["DimuonMassfinal"])
    rdf = rdf.Define("Quadruplet_Mass_eq","BsJPsiPhiMass(Dimu_OS_max, Dimu_OS_min, RefittedSV_Mass)")
    rdf = rdf.Define("isJPsiPhi","BsJPsiPhi(Dimu_OS_max, Dimu_OS_min)")
    return rdf

def DiMuVar_2(rdf, branches):
    rdf = rdf.Define("Dimu_combinations","Dimuon_v2(mu_index, Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt, MuonPt, MuonCharge)")
    rdf = rdf.Define("OS_mass","Vtx_quantity(Dimu_combinations, Vtx12_mass, Vtx23_mass, Vtx13_mass, Vtx14_mass, Vtx24_mass, Vtx34_mass)")
    rdf = rdf.Define("OS_mass_err","Vtx_quantity(Dimu_combinations, Vtx12_mass_err, Vtx23_mass_err, Vtx13_mass_err, Vtx14_mass_err, Vtx24_mass_err, Vtx34_mass_err)")
    rdf = rdf.Define("OS_Chi2","Vtx_quantity(Dimu_combinations, Vtx12_Chi2, Vtx23_Chi2, Vtx13_Chi2, Vtx14_Chi2, Vtx24_Chi2, Vtx34_Chi2)")

    for i in range(1, 3):  # i -> 1 to 2
        for j in range(1, 3):  # j -> 1 to 2
            branch_name = f"OS{i}v{j}"
            branches.append(branch_name+"_mass")
            rdf = rdf.Define(branch_name+"_mass", f"flattering(OS_mass, {(i-1) + 2*(j-1)})")
            branches.append(branch_name+"_mass_err")
            rdf = rdf.Define(branch_name+"_mass_err", f"flattering(OS_mass_err, {(i-1) + 2*(j-1)})")
            branches.append(branch_name+"_Chi2")
            rdf = rdf.Define(branch_name+"_Chi2", f"flattering(OS_Chi2, {(i-1) + 2*(j-1)})")
    
    branches.append("Dimu_OS_max")
    branches.append("Dimu_OS_min")
    branches.append("Quadruplet_Mass_eq")
    branches.append("Quadruplet_Mass_eq_old")
    branches.append("isJPsiPhi")
    rdf = rdf.Define("DimuonMassfinal","DimuonMassfinal(OS1v1_mass, OS2v1_mass, OS1v2_mass, OS2v2_mass)")
    rdf = rdf.Define("Dimu_OS_max", flat0D_double(0), ["DimuonMassfinal"])
    rdf = rdf.Define("Dimu_OS_min", flat0D_double(1), ["DimuonMassfinal"])
    rdf = rdf.Define("Quadruplet_Mass_eq_old","BsJPsiPhiMass(Dimu_OS_max, Dimu_OS_min, Quadruplet_Mass)")
    rdf = rdf.Define("Quadruplet_Mass_eq","BsJPsiPhiMass(Dimu_OS_max, Dimu_OS_min, RefittedSV_Mass)")
    rdf = rdf.Define("isJPsiPhi","BsJPsiPhi(Dimu_OS_max, Dimu_OS_min)")

    return rdf

def DiMassVar_control(rdf, branches, analysis_type):
    branches.append("Dimu_mass")
    branches.append("Ditrk_mass")
    rdf = rdf.Define("Di_mass", "DiMass"+analysis_type+"(Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt, Mu3_Eta, Mu4_Eta, Mu3_Phi, Mu4_Phi, MuonPt, MuonEta, MuonPhi, MuonEnergy)")
    rdf = rdf.Define("Dimu_mass", flat0D_double(0), ["Di_mass"])
    rdf = rdf.Define("Ditrk_mass", flat0D_double(1), ["Di_mass"])
    return rdf

    
def Gen_ct(rdf, branches, analysis_type, isMC):
    branches.append("Gen_ct_control")
    branches.append("Gen_ct_signal")
    if(isMC>0):
        if(analysis_type=="B4mu"):
            rdf = rdf.Define("Gen_ct_control", "Gen_ct(\"contol4mu\" ,MuonPt, MuonEta, MuonPhi, Mu1_Pt, Mu1_Eta, Mu1_Phi, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi, GenParticle_Pt, GenParticle_Eta, GenParticle_Phi, GenParticle_PdgId, GenParticle_MotherPdgId, GenParticle_GrandMotherPdgId, GenParticle_vx, GenParticle_vy, GenParticle_vz)")
            rdf = rdf.Define("Gen_ct_signal", "Gen_ct(\"signal\" ,MuonPt, MuonEta, MuonPhi, Mu1_Pt, Mu1_Eta, Mu1_Phi, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi, GenParticle_Pt, GenParticle_Eta, GenParticle_Phi, GenParticle_PdgId, GenParticle_MotherPdgId, GenParticle_GrandMotherPdgId, GenParticle_vx, GenParticle_vy, GenParticle_vz)")
        if(analysis_type=="B2mu2K"):
            rdf = rdf.Define("Gen_ct_signal", "1.0")
            rdf = rdf.Define("Gen_ct_control", "Gen_ct(\"contol2mu\" ,MuonPt, MuonEta, MuonPhi, Mu1_Pt, Mu1_Eta, Mu1_Phi, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi, GenParticle_Pt, GenParticle_Eta, GenParticle_Phi, GenParticle_PdgId, GenParticle_MotherPdgId, GenParticle_GrandMotherPdgId, GenParticle_vx, GenParticle_vy, GenParticle_vz)")
        if(analysis_type=="B2muKpi"):
            rdf = rdf.Define("Gen_ct_signal", "1.0")
            rdf = rdf.Define("Gen_ct_control", "1.0")
    if(isMC==0):
        rdf = rdf.Define("Gen_ct_control", "1.0")
        rdf = rdf.Define("Gen_ct_signal", "1.0")

    return rdf
    
def GenVar(rdf, branches, isMC):
    if isMC != 0:
        rdf = rdf.Define("gen_info", "GenMatching_v2(mu_index, MuonPt, MuonEta, MuonPhi, Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt, GenParticle_Pt, GenParticle_Pt_v2, GenParticle_Eta_v2, GenParticle_Phi_v2,  GenParticle_PdgId, GenParticle_MotherPdgId, GenParticle_GrandMotherPdgId)")
        rdf = rdf.Define("GenMu_Pt", flat1D_double(0), ["gen_info"])
        rdf = rdf.Define("GenMu_Eta", flat1D_double(1), ["gen_info"])
        rdf = rdf.Define("GenMu_Phi", flat1D_double(2), ["gen_info"])
    for mu in range(1,5):
        branches.append("GenMu"+str(mu)+"_Pt")
        branches.append("GenMu"+str(mu)+"_Eta")
        branches.append("GenMu"+str(mu)+"_Phi")
        if isMC != 0:
            rdf = rdf.Define("GenMu"+str(mu)+"_Pt", flat0D_double(mu-1), ["GenMu_Pt"])
            rdf = rdf.Define("GenMu"+str(mu)+"_Eta", flat0D_double(mu-1), ["GenMu_Eta"])
            rdf = rdf.Define("GenMu"+str(mu)+"_Phi", flat0D_double(mu-1), ["GenMu_Phi"])
        else:
            rdf = rdf.Define("GenMu"+str(mu)+"_Pt", add_double(0.))
            rdf = rdf.Define("GenMu"+str(mu)+"_Eta", add_double(0.))
            rdf = rdf.Define("GenMu"+str(mu)+"_Phi", add_double(0.))
    return rdf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="--Analisis inputs: directory, number of files, index")
    parser.add_argument("--index", type=int, help="index for condor submission")
    parser.add_argument("--delta", type=int, help="Number of files per submission")
    parser.add_argument("--directory_IN", type=str, help="Root files directory")
    parser.add_argument("--directory_OUT", type=str, help="Output directory")
    parser.add_argument("--isMC", type=int, help="0 for data 1 for MC")
    parser.add_argument("--analysis_type", type=str, help="Analysis type")
    args = parser.parse_args()
    index = args.index
    delta = args.delta
    directory = args.directory_IN
    output_dir = args.directory_OUT
    isMC = args.isMC
    analysis_type = args.analysis_type
    if analysis_type=="Norm":
        analysis_type="B4mu"

    if not output_dir.endswith("/"):
        output_dir= output_dir + "/"

    print(time.ctime(time.time()), " -- Starting!")
    file_root = list_of_root_files(directory)
    selected_files = select_root_files(file_root, index , delta)
    if not selected_files: print("The vector is empty. End execution."); exit()

    tree_dir_name = {
        "B2mu2K": "TreeB2mu2K",
        "B2muKpi": "TreeB2muKpi",
        "B4mu": "TreeMakerBkg"
    }.get(analysis_type) or (print("Wrong analysis type. End execution.") or exit())
    
    start_2 = time.time()
    rdf = RDataFrame(tree_dir_name+"/ntuple", selected_files) # Load data
    print(time.ctime(time.time()), " -- Data loaded!")
    
    #Find best Quadruplet
    rdf = rdf.Define("isMC", add_int(isMC))
    branches=["evt", "isMC", "run", "lumi", "nPileUpInt", "PVCollection_Size"]
    rdf.Snapshot("FinalTree", output_dir + "Analyzed_Data_index_"+str(index)+".root", branches)
    exit()

    if(analysis_type=="B4mu"):
        rdf = rdf.Define("Quadruplet_indexs","B4mu_QuadSel(isMC, evt, MuonPt, MuonEta, MuonPhi, RefTrack1_Pt, Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt, NGoodQuadruplets, QuadrupletVtx_Chi2, RefittedSV_Mass, Muon_isGlobal, Muon_isPF, Muon_isLoose, Muon_isMedium, Muon_isTight, Muon_isSoft, MuonPt_HLT, MuonEta_HLT, MuonPhi_HLT, FlightDistBS_SV_Significance, Muon_vz, GenParticle_Pt, GenParticle_Eta, GenParticle_Phi, GenParticle_Pt_v2, GenParticle_Eta_v2, GenParticle_Phi_v2, GenParticle_PdgId, GenParticle_MotherPdgId, GenParticle_GrandMotherPdgId, vtx_prob, QuadrupletVtx_x, QuadrupletVtx_y, RefittedPV_x, RefittedPV_y, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi, Quadruplet_Charge)")
    else:
        rdf = rdf.Define("remove_duplicate",analysis_type+"_CombSel(Mu3_Pt, Mu4_Pt, Mu3_Eta, Mu4_Eta, Mu3_Phi, Mu4_Phi, QuadrupletVtx_Chi2)")
        rdf = rdf.Define("Quadruplet_indexs","B2muX_QuadSel(remove_duplicate, isMC, evt, MuonPt, MuonEta, MuonPhi, RefTrack1_Pt, Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt, Mu3_Eta, Mu4_Eta, NGoodQuadruplets, QuadrupletVtx_Chi2, RefittedSV_Mass, Muon_isGlobal, Muon_isPF, Muon_isLoose, Muon_isMedium, Muon_isTight, Muon_isSoft, MuonPt_HLT, MuonEta_HLT, MuonPhi_HLT, FlightDistBS_SV_Significance, Muon_vz, GenParticle_Pt, GenParticle_Pt_v2, GenParticle_Eta_v2, GenParticle_Phi_v2, GenParticle_PdgId, GenParticle_MotherPdgId, GenParticle_GrandMotherPdgId, vtx_prob, QuadrupletVtx_x, QuadrupletVtx_y, RefittedPV_x, RefittedPV_y, Quadruplet_Pt, Quadruplet_Eta, Quadruplet_Phi)")
    
    rdf = rdf.Define("Quadruplet_index", flat0D_int(0), ["Quadruplet_indexs"])
    rdf = rdf.Filter("Quadruplet_index>-1")
    
    rdf = Flat_MuVar(rdf, branches) #Flat muon pt eta phi
    if(analysis_type=="B4mu"):
        rdf = rdf.Define("mu_index", "get_4index(MuonPt, Mu1_Pt, Mu2_Pt, Mu3_Pt, Mu4_Pt)")
        rdf = rdf.Define("dz_max", "DeltaZmax(mu_index, Muon_vz)")
        branches.append("dz_max")
    else:
        rdf = rdf.Define("mu_index", "get_2index(MuonPt, Mu1_Pt, Mu2_Pt)")
    
    if(analysis_type=="B4mu"):
        rdf, branches = MuonIDs(rdf, branches) #Add muonIDs
    else:
        rdf, branches = MuonIDs(rdf, branches, n_muons=2) #Add muonIDs

    rdf, vertex_chi2 = QuadMuVar(rdf, branches, analysis_type) #Quadruplet variables
    rdf = MVA_inputs(rdf, branches) #Define MVA input variables
    if(analysis_type=="B4mu"):
        #rdf = DiMuVar(rdf, branches, vertex_chi2) #Define Di-Muon variables
        rdf = DiMuVar_2(rdf, branches) #Define Di-Muon variables
        #rdf, branches = HLT_quantities(rdf, branches)
        #rdf = Gen_ct(rdf, branches, analysis_type, isMC)
        #rdf = GenVar(rdf, branches, isMC) #Gen-Level variables for control channel
    
    rdf.Snapshot("FinalTree", output_dir + "Analyzed_Data_index_"+str(index)+".root", branches)
    
    print(time.ctime(time.time()), " -- Performed ",rdf.GetNRuns()," loops")
    end = time.time()

    del rdf
    del branches

    print('Partial execution time ', end-start_2)
    print('Total execution time ', end-start)
