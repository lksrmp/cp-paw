import argparse
from argparse import RawTextHelpFormatter
import sys
import numpy as np
import matplotlib.pyplot as plt

DPI = 150  # dots per inch
TRANSPARENT = True  # make background transparent
FIG_SIZE = (6.4, 4.8)  # figure size in inches
FONTSIZE = 14  # font size
XLIM = None  # x-axis limits (min, max)
YLIM = None  # y-axis limits (min, max)
XTICKS = None  # y-axis ticks [tick1, tick2, ...]
XTICK_LABELS = None  # y-axis tick labels e.g. [r"$\Gamma$", r"$X$", r"$W$", ...]
COLORS = None  # colors for the bands e.g. ['black']*4 + ['blue']*4


def plot_data(kgrid, bands, filename):
    """
    Plot the data using matplotlib.

    Parameters:
    kgrid (array-like): The x-axis values.
    bands (array-like): The y-axis values.
    filename (str): The filename to save the plot. If None, the plot will be displayed.

    Returns:
    None
    """
    plt.rcParams.update({"font.size": FONTSIZE})
    fig, ax = plt.subplots(layout="constrained", dpi=DPI)
    # ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useOffset=False)
    lines = []
    if COLORS:
        if len(COLORS) != bands.shape[1]:
            print("COLORS and bands must have the same length.")
            print(f"COLORS: {len(COLORS)}, bands: {bands.shape[1]}")
            sys.exit(1)
        colorlist = COLORS
    else:
        colorlist = ['black']*bands.shape[1]
    for i in range(bands.shape[1]):
        (line,) = ax.plot(kgrid, bands[:, i], color=colorlist[i])
        lines.append(line)
    ax.grid()
    ax.set_ylabel('energy [eV]')
    if XTICKS:
        ax.set_xticks(XTICKS)
    else:
        ax.set_xticks([])
    if XTICK_LABELS and XTICKS:
        if(len(XTICK_LABELS) != len(XTICKS)):
            print("XTICKS and XTICK_LABELS must have the same length.")
            sys.exit(1)
        ax.set_xticklabels(XTICK_LABELS)
    if XLIM:
        ax.set_xlim(XLIM)
    else:
        ax.set_xlim(np.min(kgrid[0]), np.max(kgrid[-1]))
    if YLIM:
        ax.set_ylim(YLIM)

    if filename:
        if filename.endswith(".png"):
            fig.savefig(filename, format="png", transparent=TRANSPARENT)
        elif filename.endswith(".jpg"):
            fig.savefig(filename, format="jpg")
        else:
            print("Unknown file extension. Allowed are .eps and .jpg.")
            sys.exit(1)
    else:
        plt.show()


def parse_args():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Plot time-dependent data from protocol file.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        help=(
            "Print to file rather than starting xmgrace.\n"
            "Extensions eps,png produce graphics files,\n"
            "others or no extensions produce xmgrace input files.\n"
            "Process with xmgrace -nxy filename."
        ),
    )
    # parser.add_argument('-b', '--batch', type=str,
    #                     help="Attach batch file to xmgrace settings.")
    parser.add_argument(
        "file", type=str, help="file containing band data"
    )
    return parser.parse_args()

def read_values_to_array(filepath):
    """
    Read values from a file and return them as a 2D array.

    Parameters:
    filepath (str): The path to the file.

    Returns:
    numpy.ndarray or None: The 2D array containing the values read from the file, or None if an error occurred.

    Raises:
    OSError: If the file is not found.
    ValueError: If the file could not be read.

    """
    try:
        data = np.loadtxt(filepath)
    except OSError:
        print(f"File {filepath} not found.")
        return None
    except ValueError:
        print(f"File {filepath} could not be read.")
        return None
    if data.ndim != 2:
        print(f"File {filepath} has wrong dimension.")
        return None
    if data.shape[0] < 2:
        print(f"File {filepath} has too few rows.")
        return None
    if data.shape[1] < 2:
        print(f"File {filepath} has too few columns.")
        return None
    return data


def main():
    """
    Main function that reads values from a file and plots the data.

    Args:
        None

    Returns:
        None
    """
    args = parse_args()
    
    # Read the values from the file
    values = read_values_to_array(args.file)
    if values is None:
        sys.exit(1)
    kgrid = values[:, 0]
    bands = values[:, 1:]
    plot_data(kgrid, bands, args.out)


if __name__ == "__main__":
    main()