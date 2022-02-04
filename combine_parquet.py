import pandas as pd
import json
import argparse
import os

def loadSummaries(args):
  summaries = []
  for summary_path in args.summary_input:
    with open(summary_path, "r") as f:
      summaries.append(json.load(f)["sample_id_map"])

  return summaries

def loadDataFrames(args):
  dfs = []
  for parquet_path in args.parquet_input:
    dfs.append(pd.read_parquet(parquet_path))

  return dfs

def invertSummary(summary):
  inverted_summary = {}
  for proc_name in summary.keys():
    inverted_summary[summary[proc_name]] = proc_name
  return inverted_summary

def idToName(dfs, summaries):
  for i, df in enumerate(dfs):
    inverted_summary = invertSummary(summaries[i])
    df.process_id.replace(inverted_summary, inplace=True)  
  return dfs

def nameToId(df, summary):
  df.process_id.replace(summary, inplace=True)
  return df

def removeExcludedProcesses(dfs, summaries, args):
  for i, summary in enumerate(summaries):
    for process in list(summary.keys()):
      if process in args.exclude_procs:
        dfs[i] = dfs[i][dfs[i].process_id != summary[process]]
        del summary[process]        

  return dfs, summaries

def mergeDataFrames(args):
  summaries = loadSummaries(args)
  dfs = loadDataFrames(args)

  dfs, summaries = removeExcludedProcesses(dfs, summaries, args)
  dfs = idToName(dfs, summaries)

  merged_summary = {}
  for summary in summaries: merged_summary.update(summary)

  if len(merged_summary.keys()) == sum([len(summary.keys()) for summary in summaries]): #no conflict
    print(">> No overlapping process names")
    merged_df = pd.concat(dfs)
    merged_summary = {process:i for i, process in enumerate(merged_df.process_id.unique())}
  else:
    print("Have not implemented a way to deal with conflicting process names yet. Script will now exit.")
    exit()

  merged_df = nameToId(merged_df, merged_summary)

  merged_df.to_parquet(args.parquet_output)
  with open(args.summary_output, "w") as f:
    json.dump({"sample_id_map": merged_summary}, f, indent=4)  

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--parquet-input', '-p', type=str, nargs='+', help="List of parquet files to merge, e.g. low_mass.parquet sm_trigger.parquet ...")
  parser.add_argument('--summary-input', '-s', type=str, nargs='+', help="List of summary json files to merge, e.g. low_mass_summary.json sm_trigger_summary.json ...")

  parser.add_argument('--parquet-output', type=str, default="merged_nominal.parquet")
  parser.add_argument('--summary-output', type=str, default="summary.json")

  parser.add_argument('--exclude-procs', '-e', type=str, nargs='+', default=[], help="List of processes to not include in the merging, e.g. Diphoton TTGamma ...")

  parser.add_argument('--force', '-f', default=False, action="store_true", help="Overwrite output parquet and summary files without asking permission.")

  args = parser.parse_args()

  if os.path.exists(args.parquet_output):
    if input("%s already exists. Should we continue anyway? (y/n): "%args.parquet_output) == "n":
      exit()
  if os.path.exists(args.summary_output):
    if input("%s already exists. Should we continue anyway? (y/n): "%args.summary_output) == "n":
      exit()

  assert len(args.parquet_input) == len(args.summary_input)

  mergeDataFrames(args)

