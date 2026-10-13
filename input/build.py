import json
from subprocess import run

sample_name = "DYtoMuTauh_MLL_50_madgraphMLM"
sample = "/DYto2Tau-4Jets_Bin-MLL-50_Fil-MuTauh_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM"
command = f"/cvmfs/cms.cern.ch/common/dasgoclient --query='file dataset={sample}'"

file_list = run(command, shell=True, capture_output=True, text=True).stdout.splitlines()
prefix = "root://cms-xrd-global.cern.ch/"
final_file_list = [f"{prefix}{file}" for file in file_list]

output = {
    sample_name: final_file_list
}

with open("input/2024_DYtoMuTauh.json", "w") as f:
    json.dump(output, f, indent=4)

