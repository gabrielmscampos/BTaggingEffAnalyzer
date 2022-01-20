import argparse
import json
import os

import anatools.analysis as ana
import anatools.data as data
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import btaggingeffmaps as btef

ana.start()

KEYS = [
    "dataset_name",
    "max_unc",
    "eta_bins",
    "year",
    "apv",
    "algo",
    "working_point",
    "pt_min",
    "pt_max",
    "step_size",
    "unc_stop",
    "unc_increase",
    "find_best_unc",
]


def find_maximum_pt_cut(_datasets, _pt_max=1000.0, threshold=0.001):
    """
    Find maximum pt cut in dataset in which excluded events are lesser then (threshold*100)%

    Args:
        datasets (dict[pd.DataFrame]): Dict of datasets in ANATools format
        pt_max (float): Initial bet for maximum pt cut
        threshold (float): percentage of number of events accepted after cut

    Returns:
        float: Maximum pt selection criteria
    """
    stop_search = False
    while stop_search is False:
        above_thr = False
        for _dataset_name, _dataset in _datasets.items():
            df_alljets = _dataset.reset_index(drop=True)
            n_total = df_alljets.evtWeight.sum()
            n_above_cut = df_alljets[df_alljets.Jet_pt > _pt_max].evtWeight.sum()
            quant = n_above_cut / n_total
            if quant > threshold:
                above_thr = True
                quant = "{:.1f}%".format(100 * quant)
                print(_dataset_name, "->", quant)
        if above_thr:
            _pt_max += 10
        else:
            stop_search = True
    print(f"Chosen pt_max = {_pt_max}")
    return _pt_max


if __name__ == "__main__":

    #######################
    # Argument parser
    #######################
    parser = argparse.ArgumentParser()
    parser.add_argument("--period", dest="period", type=str)
    parser.add_argument("--basedir", dest="basedir", type=str)
    parser.add_argument("--year", dest="year", type=str)
    parser.add_argument("--apv", dest="apv", action="store_true", required=False)
    parser.add_argument("--algo", dest="algo", type=str)
    parser.add_argument("--working-point", dest="working_point", type=str)
    parser.add_argument("--eta-bins", dest="eta_bins", type=float, nargs="+")
    parser.add_argument("--pt-min", dest="pt_min", type=float)
    parser.add_argument(
        "--pt-max", dest="pt_max", type=float, required=False, default=1000.0
    )
    parser.add_argument(
        "--pt-max-thr", dest="pt_max_thr", type=float, required=False, default=0.001
    )
    parser.add_argument("--step-size", dest="step_size", type=float)
    parser.add_argument("--unc-stop", dest="unc_stop", type=float)
    parser.add_argument("--unc-increase", dest="unc_increase", type=float)
    parser.add_argument(
        "--not-find-best-unc",
        dest="find_best_unc",
        action="store_false",
        required=False,
    )
    parser.add_argument(
        "--output-path",
        dest="output_path",
        type=str,
        required=False,
        default="./output",
    )

    parser.set_defaults(apv=False)
    parser.set_defaults(find_best_unc=True)
    args = parser.parse_args()

    print(
        f"\nperiod: {args.period}\n"
        f"basedir: {args.basedir}\n"
        f"apv: {args.apv}\n"
        f"algo: {args.algo}\n"
        f"working_point: {args.working_point}\n"
        f"eta_bins: {args.eta_bins}\n"
        f"pt_min: {args.pt_min}\n"
        f"pt_max: {args.pt_max}\n"
        f"pt_max_thr: {args.pt_max_thr}\n"
        f"step_size: {args.step_size}\n"
        f"unc_stop: {args.unc_stop}\n"
        f"unc_increase: {args.unc_increase}\n"
        f"find_best_unc: {args.find_best_unc}\n"
        f"output_path: {args.output_path}"
    )

    year = "20" + args.period

    #######################
    # Read datasets
    #######################
    print("\nReading datasets...")
    datasets = data.read_files(args.basedir, args.period)

    DY = [
        datasets["DYJetsToLL_Pt-Inclusive"],
    ]

    datasets["DYJetsToLL"] = pd.concat(DY).reset_index(drop=True)

    for dtname in list(datasets.keys()):
        if dtname.startswith("DYJetsToLL_"):
            del datasets[dtname]

    #######################
    # Find maximum pt cut
    #######################
    print("\nFinding maximum pt cut...")
    pt_max = find_maximum_pt_cut(datasets, args.pt_max, args.pt_max_thr)

    #######################
    # Make efficiency maps
    #######################
    print(f"\nMaking efficiency maps ({len(datasets.keys())} datasets)...")
    efficiency_map_fname = (
        f"btageffmap-{args.algo}-{args.working_point}-{year}-{args.apv}.json"
    )
    uncertainty_map_fname = (
        f"btaguncmap-{args.algo}-{args.working_point}-{year}-{args.apv}.json"
    )
    efficiency_map_fpath = os.path.join(args.output_path, efficiency_map_fname)
    uncertainty_map_fpath = os.path.join(args.output_path, uncertainty_map_fname)

    # Use uncertainty map if already exists, otherwise create with default values
    if os.path.exists(uncertainty_map_fpath) is False:
        accepted_unc = {
            dataset_name: {"b": 0.001, "c": 0.001, "udsg": 0.001}
            for dataset_name in datasets
        }
    else:
        with open(uncertainty_map_fpath, "r", encoding="utf-8") as f:
            accepted_unc = json.load(f)

    eff_maps = {}
    unc_maps = {}
    for dataset_number, (dataset_name, dataset) in enumerate(datasets.items(), 1):
        print(f"    - ({dataset_number}/{len(datasets.keys())}) {dataset_name}")
        df = dataset.reset_index(drop=True)
        mmap = btef.BTaggingEfficiencyMap(df, args.eta_bins)
        mmap.calib(
            year,
            args.apv,
            args.algo,
            args.working_point,
        )
        eff_map, unc_map = mmap.make(
            args.pt_min,
            args.pt_max,
            args.step_size,
            accepted_unc.get(dataset_name),
            args.find_best_unc,
            args.unc_stop,
            args.unc_increase,
        )
        eff_maps[dataset_name] = eff_map
        unc_maps[dataset_name] = unc_map

    with open(efficiency_map_fpath, "w", encoding="utf-8") as f:
        json.dump(eff_maps, f)

    with open(uncertainty_map_fpath, "w", encoding="utf-8") as f:
        json.dump(unc_maps, f, indent=4, ensure_ascii=False)

    #######################
    # Plot effiency per eta bin
    #######################
    print("\nPlotting efficiency per eta bin for each dataset efficiency map...")
    outname_year = "APV_" + year if args.apv else year
    outpath_year = os.path.join(args.output_path, outname_year)

    if os.path.isdir(outpath_year) is False:
        os.makedirs(outpath_year)

    for dataset_name, eff_map in eff_maps.items():
        plot_fpath = os.path.join(
            outpath_year,
            f"{dataset_name}_effetabin_{args.algo}-{args.working_point}.png",
        )
        fig = plt.figure(figsize=(24, 8))
        for idx, hf in enumerate(eff_map.keys()):
            dt = eff_map.get(hf)
            eta_min = np.array([d.get("eta_min") for d in dt])
            pt_min = np.array([d.get("pt_min") for d in dt])
            eff = np.array([d.get("eff") for d in dt])

            shape = (int(pt_min.size / 3), 3)
            eta_min = eta_min.reshape(shape)[0]
            pt_min = pt_min.reshape(shape)[:, 0]
            eff = eff.reshape(shape)

            ax = fig.add_subplot(1, 3, idx + 1)
            for jdx, eta in enumerate(eta_min):
                ax.plot(pt_min, eff[:, jdx], label=f"eta bin {eta}")

            ax.set_title(f"{dataset_name}: {hf}-jets", fontsize=18)
            ax.set_xlabel("pt")
            ax.set_ylabel("efficiency")
            ax.legend()

        fig.tight_layout()
        plt.savefig(plot_fpath, facecolor="white", dpi=200)
        plt.close()

    print("\nRoutine finished.")
