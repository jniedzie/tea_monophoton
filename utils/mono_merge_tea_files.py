import os
import argparse

from lbl_paths import base_path, processes, skim


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", action="store_true")
    return parser.parse_args()


def main():
    args = get_args()

    merge_histograms = not args.trees

    for process in processes:
        hist_path = "histograms" if merge_histograms else ""
        input_path = f"{base_path}/{process}/{skim}/{hist_path}/*.root"

        hist_path = "_histograms" if merge_histograms else ""
        output_path = f"{base_path}/{process}/merged_{skim}{hist_path}.root"

        print(f"{process=}")
        print(f"{output_path=}")

        os.system(f"rm {output_path}")
        os.system(f"hadd -f -j -k {output_path} {input_path}")


if __name__ == "__main__":
    main()
