import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep
mplhep.set_style("CMS")
plt.rcParams["figure.figsize"] = (12.5,10)
#plt.rcParams['figure.constrained_layout.use'] = True

import pandas as pd
import numpy as np

import argparse
import os

import json
from collections import OrderedDict

from tqdm import tqdm

bkg_groups = {
  'Diphoton': ["Diphoton_MGG-40to80", "Diphoton_MGG-80toInf"],
  'GJets': ['GJets_HT-100To200', 'GJets_HT-200To400', 'GJets_HT-400To600', 'GJets_HT-40To100', 'GJets_HT-600ToInf'],
  'TT': ['TTGG', 'TTGamma', 'TTJets'],
  'SM Higgs': ['VBFH_M125', 'VH_M125', 'ggH_M125', 'ttH_M125'],
  'VGamma': ['WGamma', 'ZGamma']
}
bkg_processes = []
for group in bkg_groups.keys(): bkg_processes.extend(bkg_groups[group])

columns_to_plot = ["Diphoton_mass", "Diphoton_pt", "Diphoton_eta", "Diphoton_phi", "LeadPhoton_pt", "LeadPhoton_eta", "LeadPhoton_phi", "LeadPhoton_mass", "LeadPhoton_mvaID", "LeadPhoton_genPartFlav", "LeadPhoton_pixelSeed", "SubleadPhoton_pt", "SubleadPhoton_eta", "SubleadPhoton_phi", "SubleadPhoton_mass", "SubleadPhoton_mvaID", "SubleadPhoton_genPartFlav", "SubleadPhoton_pixelSeed", "jet_1_pt", "jet_1_eta", "jet_1_phi", "jet_1_mass", "jet_1_btagDeepFlavB", "jet_2_pt", "jet_2_eta", "jet_2_phi", "jet_2_mass", "jet_2_btagDeepFlavB", "b_jet_1_btagDeepFlavB", "b_jet_2_btagDeepFlavB", "tau_candidate_1_pt", "tau_candidate_1_eta", "tau_candidate_1_phi", "tau_candidate_1_mass", "tau_candidate_1_charge", "tau_candidate_1_id", "tau_candidate_2_pt", "tau_candidate_2_eta", "tau_candidate_2_phi", "tau_candidate_2_mass", "tau_candidate_2_charge", "tau_candidate_2_id", "tau_candidate_3_pt", "tau_candidate_3_eta", "tau_candidate_3_phi", "tau_candidate_3_mass", "tau_candidate_3_charge", "tau_candidate_3_id", "n_jets", "n_leptons", "n_electrons", "n_muons", "n_taus", "n_iso_tracks", "ditau_pt", "ditau_eta", "ditau_phi", "ditau_mass", "ditau_dR", "ditau_lead_lepton_pt", "ditau_lead_lepton_eta", "ditau_lead_lepton_phi", "ditau_lead_lepton_mass", "ditau_lead_lepton_id", "ditau_lead_lepton_charge", "ditau_sublead_lepton_pt", "ditau_sublead_lepton_eta", "ditau_sublead_lepton_phi", "ditau_sublead_lepton_mass", "ditau_sublead_lepton_id", "ditau_sublead_lepton_charge", "dilep_leadpho_mass", "dilep_subleadpho_mass"]

"""
#pastel
colour_schemes = {
  4: ['#b3e2cd','#fdcdac','#cbd5e8','#f4cae4'],
  5: ['#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9'],
  6: ['#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae'],
  7: ['#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae','#f1e2cc'],
  8: ['#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae','#f1e2cc','#cccccc']
}
"""
colour_schemes = {
  4: ['#a6cee3','#1f78b4','#b2df8a','#33a02c'],
  5: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99'],
  6: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c'],
  7: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f'],
  8: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00']
}

def preprocess(df):
  df = dataFrameBoolToInt(df)

  pixel_veto = (df.LeadPhoton_pixelSeed==0) & (df.SubleadPhoton_pixelSeed==0)
  df = df[pixel_veto]
  return df

def dataFrameBoolToInt(df):
  type_conversion = {}
  for column in df.columns:
    if df[column].dtype == bool:
      type_conversion[column] = int
  return df.astype(type_conversion)

def createDefaultConfig(data, bkg, sig):
  config = OrderedDict()
  for column in columns_to_plot:
    d = data[column][data[column]!=-999]
    b = bkg[column][bkg[column]!=-999]
    s = sig[column][sig[column]!=-999]

    #low = min([d.min(), b.min(), s.min()])
    #high = max([d.max(), b.max(), s.max()])
    low = min([d.quantile(0.05), b.quantile(0.05), s.quantile(0.05)])
    high = max([d.quantile(0.95), b.quantile(0.95), s.quantile(0.95)])
    config[column] = {"range": [float(low),float(high)]}
  return config

def writeDefaultConfig(data, bkg, sig):
  config = createDefaultConfig(data, bkg, sig)
  with open("config.json", "w") as f:
    json.dump(config, f, indent=4)

def getConfig(data, bkg, sig, args):
  if args.config != None:
    with open(args.config, "r") as f:
      cfg = json.load(f)
  else:
    cfg = createDefaultConfig(data, bkg, sig)
    
  return cfg

def createBkgStack(bkg, column, group=True):
  bkg_stack = []
  bkg_stack_w = []
  bkg_stack_labels = []

  if not group:
    for proc in bkg_processes:
      bkg_stack.append(bkg[bkg.process_id==proc_dict[proc]][column])
      bkg_stack_w.append(bkg[bkg.process_id==proc_dict[proc]]["weight_central"])
      bkg_stack_labels.append(proc)
  else:
    for bkg_group in bkg_groups.keys():
      proc_ids = [proc_dict[proc] for proc in bkg_groups[bkg_group]]
      bkg_stack.append(bkg[bkg.process_id.isin(proc_ids)][column])
      bkg_stack_w.append(bkg[bkg.process_id.isin(proc_ids)]["weight_central"])
      bkg_stack_labels.append(bkg_group)

  is_sorted = False
  while not is_sorted:
    is_sorted = True
    for i in range(len(bkg_stack)-1):
      if bkg_stack_w[i].sum() > bkg_stack_w[i+1].sum():
        is_sorted = False
        bkg_stack[i], bkg_stack[i+1] = bkg_stack[i+1], bkg_stack[i]
        bkg_stack_w[i], bkg_stack_w[i+1] = bkg_stack_w[i+1], bkg_stack_w[i]
        bkg_stack_labels[i], bkg_stack_labels[i+1] = bkg_stack_labels[i+1], bkg_stack_labels[i]

  return bkg_stack, bkg_stack_w, bkg_stack_labels

def decayToMath(channel):
  if channel == "gg":
    return r"\gamma\gamma"
  else:
    return r"\tau\tau"

def getSigLabel(args):
  if "NMSSM" in args.sig_proc:
    split_name = args.sig_proc.split("_")
    Y_decay = decayToMath(split_name[3])
    H_decay = decayToMath(split_name[5])
    X_mass = int(split_name[7])
    Y_mass = int(split_name[9])
    label = r"$X_{%d} \rightarrow Y_{%d}(\rightarrow %s)  H(\rightarrow %s)$"%(X_mass, Y_mass, Y_decay, H_decay)
  elif "radion" in args.sig_proc:
    X_mass = int(args.sig_proc.split("M")[1].split("_")[0])
    label = r"$X_{%d} \rightarrow HH \rightarrow \gamma\gamma\tau\tau$"%X_mass
  else:
    label = args.sig_proc
  return label

def adjustLimits(x, ys, ax):
  data_to_display = ax.transData.transform
  display_to_data = ax.transData.inverted().transform

  tx = lambda x: data_to_display((x,0))[0]
  tx_inv = lambda x: display_to_data((x,0))[0]
  ty = lambda x: data_to_display((0,x))[1]
  ty_inv = lambda x: display_to_data((0,x))[1]

  xlow, xhigh = tx(ax.get_xlim()[0]), tx(ax.get_xlim()[1])
  ylow, yhigh = ty(ax.get_ylim()[0]), ty(ax.get_ylim()[1])
  
  #top side
  ybound = ylow + (yhigh-ylow)*0.60
  max_y = np.array(ys).max()
  top_distance_to_move = ty(max_y) - ybound
  
  #right side
  xbound = xlow + (xhigh-xlow)*0.75
  ybound = ylow + (yhigh-ylow)*0.20
  max_y = np.array(ys).T[x>tx_inv(xbound)].max()
  right_distance_to_move = ty(max_y) - ybound

  if right_distance_to_move <= 0:
    ax.legend(ncol=1, loc="upper right", markerfirst=False)
  elif right_distance_to_move < top_distance_to_move:
    ax.legend(ncol=1, loc="upper right", markerfirst=False)
    ax.set_ylim(top = ty_inv(yhigh + right_distance_to_move))
  else:
    ax.legend(ncol=3, loc="upper right", markerfirst=False)
    ax.set_ylim(top = ty_inv(yhigh + top_distance_to_move))

def plot(data, bkg, sig, args):  
  cfg = getConfig(data, bkg, sig, args)

  f, axs = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
  
  for column in tqdm(cfg.keys()):    
    data_hist, edges = np.histogram(data[column], bins=50, range=cfg[column]["range"], weights=data["weight_central"])
    bkg_hist, edges = np.histogram(bkg[column], bins=50, range=cfg[column]["range"], weights=bkg["weight_central"])
    sig_hist, edges = np.histogram(sig[column], bins=50, range=cfg[column]["range"], weights=sig["weight_central"])
    bin_centres = (edges[:-1]+edges[1:])/2

    bkg_stack, bkg_stack_w, bkg_stack_labels = createBkgStack(bkg, column)

    ratio = data_hist / bkg_hist
    ratio_err = np.sqrt(data_hist) / bkg_hist

    sig_sf = data_hist.sum() / sig_hist.sum()

    axs[0].hist(edges[:-1], edges, weights=sig_hist*sig_sf, label=getSigLabel(args), histtype='step', color='r', lw=3, zorder=10) #signal
    axs[0].hist(bkg_stack, edges, weights=bkg_stack_w, label=bkg_stack_labels, stacked=True, color=colour_schemes[len(bkg_stack)]) #background
    axs[0].errorbar(bin_centres, data_hist, np.sqrt(data_hist), label="Data", fmt='ko') #data
    axs[0].set_ylabel("Events")
  
    axs[1].errorbar(bin_centres, ratio, ratio_err, label="Data", fmt='ko')
    axs[1].set_xlabel(column)
    axs[1].set_ylabel("Data / MC")

    plt.sca(axs[0])
    mplhep.cms.label(llabel="Work in Progress", data=True, lumi=59, loc=0)

    adjustLimits(bin_centres, [sig_hist*sig_sf, data_hist], axs[0])
    plt.savefig(os.path.join(args.output, "%s.png"%column))
    #plt.savefig(os.path.join(args.output, "%s.pdf"%column))

    axs[0].set_yscale("log")
    axs[0].relim()
    axs[0].autoscale()
    axs[0].get_ylim()
    adjustLimits(bin_centres, [sig_hist*sig_sf, data_hist], axs[0])
    plt.savefig(os.path.join(args.output, "%s_log.png"%column))
    #plt.savefig(os.path.join(args.output, "%s_log.pdf"%column))

    axs[0].cla()
    axs[1].cla()

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--input', '-i', type=str)
  parser.add_argument('--sig-input', '-si', type=str)
  parser.add_argument('--summary', '-s', type=str)
  parser.add_argument('--sig-summary', '-ss', type=str)
  parser.add_argument('--sig-proc', '-p', type=str)
  parser.add_argument('--output', '-o', type=str, default="plots")
  parser.add_argument('--config', '-c', type=str)
  parser.add_argument('--norm', default=False, action="store_true")
  args = parser.parse_args()

  if args.sig_summary == None: args.sig_summary = args.summary

  with open(args.summary, "r") as f:
    proc_dict = json.load(f)["sample_id_map"]
  with open(args.sig_summary, "r") as f:
    sig_proc_dict = json.load(f)["sample_id_map"]

  print(">> Loading dataframes")  
  df = pd.read_parquet(args.input)
  if args.sig_input == None:
    sig_df = df
  else:
    sig_df = pd.read_parquet(args.sig_input)

  print(">> Preprocessing")
  df = preprocess(df)
  sig_df = preprocess(df)  

  data = df[df.process_id==proc_dict["Data"]]
  bkg_proc_ids = [proc_dict[bkg_proc] for bkg_proc in bkg_processes]
  bkg = df[df.process_id.isin(bkg_proc_ids)]
  sig = sig_df[sig_df.process_id==sig_proc_dict[args.sig_proc]]

  #blinding
  #data = data[(data.Diphoton_mass < 120) | (data.Diphoton_mass>130)]
  #bkg = bkg[(bkg.Diphoton_mass < 120) | (bkg.Diphoton_mass>130)]

  #normalise bkg mc to data
  if args.norm:
    bkg.loc[:, "weight_central"] = bkg.weight_central * (data.weight_central.sum() / bkg.weight_central.sum())

  os.makedirs(args.output, exist_ok=True)
   
  np.seterr(all='ignore')
  plot(data, bkg, sig, args)
  
  