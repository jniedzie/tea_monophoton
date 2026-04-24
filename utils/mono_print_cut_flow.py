import ROOT
from collections import OrderedDict

from lbl_params import luminosity, crossSections, nGenEvents, get_scale_factor
from lbl_helpers import get_cep_scale
from lbl_paths import base_path, processes, skim


def main():
    ROOT.gROOT.SetBatch(True)

    photon_sf, _ = get_scale_factor(photon=True)
    print(f"{photon_sf=}")
    
    trigger_sf = 1.0089 * 0.8716
    recoID_sf = 0.9758 * 0.9758 * 0.946 * 0.946
    che_sf = 0.9252
    nee_sf = 0.8487
    
    sf = photon_sf
    # sf = 1
    # sf = trigger_sf
    # sf = trigger_sf * recoID_sf
    # sf = trigger_sf * recoID_sf * che_sf
    # sf = trigger_sf * recoID_sf * che_sf * nee_sf
    
    for process in processes:
        input_path = f"{base_path}/{process}/merged_{skim}.root"
        print(f"Analyzing file: {input_path}")

        if process == "cep":
            scale, _ = get_cep_scale(skim)
        elif process == "collisionData" or process == "emptyBX" or process == "zeroBias":
            scale = 1
        else:
            scale = luminosity*crossSections[process]*sf
            scale /= nGenEvents[process]

        file = ROOT.TFile(input_path, "READ")
        dir = file.Get("CutFlow")
        keys = dir.GetListOfKeys()
        hist_dict = OrderedDict()

        for key in keys:
            hist = dir.Get(key.GetName())
            hist_dict[key.GetName()] = hist.GetBinContent(1)

        hist_dict = OrderedDict(
            sorted(hist_dict.items(), key=lambda x: int(x[0].split("_")[0])))

        print("CutFlow:")
        for key, value in hist_dict.items():
            print(f"{key:30}: {value:10}\t\t{value*scale:.4f}")


if __name__ == "__main__":
    main()
