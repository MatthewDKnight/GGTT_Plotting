#!/usr/bin/env bash

mkdir -p plots

for proc in NMSSM_XYH_Y_gg_H_tautau_MX_300_MY_70 NMSSM_XYH_Y_gg_H_tautau_MX_300_MY_100 NMSSM_XYH_Y_gg_H_tautau_MX_500_MY_70 NMSSM_XYH_Y_gg_H_tautau_MX_500_MY_100 NMSSM_XYH_Y_gg_H_tautau_MX_700_MY_70 NMSSM_XYH_Y_gg_H_tautau_MX_700_MY_100 NMSSM_XYH_Y_tautau_H_gg_MX_300_MY_50 NMSSM_XYH_Y_tautau_H_gg_MX_300_MY_70 NMSSM_XYH_Y_tautau_H_gg_MX_300_MY_100 NMSSM_XYH_Y_tautau_H_gg_MX_500_MY_50 NMSSM_XYH_Y_tautau_H_gg_MX_500_MY_70 NMSSM_XYH_Y_tautau_H_gg_MX_500_MY_100; do
  python -i plot.py --input low_mass.parquet --summary low_mass.json --sig-proc ${proc} -o plots/low_mass_${proc}
  python -i plot.py --input low_mass.parquet --summary low_mass.json --sig-proc ${proc} -o plots/low_mass_${proc}_norm --norm
done

for proc in radionM300_HHggTauTau radionM400_HHggTauTau radionM500_HHggTauTau radionM800_HHggTauTau radionM1000_HHggTauTau; do
  python -i plot.py --input sm_trigger.parquet --summary sm_trigger.json --sig-proc ${proc} -o plots/low_mass_${proc}
  python -i plot.py --input sm_trigger.parquet --summary sm_trigger.json --sig-proc ${proc} -o plots/low_mass_${proc}_norm --norm
done