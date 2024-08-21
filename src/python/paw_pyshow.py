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
YZOOM = None  # zoom in y-axis, only allowed if YLIM is None and one curve is plotted


def plot_data(time, plt_values, labels, xlabel, filename):
    """
    Plot data using matplotlib.

    Parameters:
    - time (array-like): The time values for the x-axis.
    - plt_values (list of array-like): The values to be plotted on the y-axis.
    - labels (list of str): The labels for each curve.
    - xlabel (str): The label for the x-axis.
    - filename (str): The filename to save the plot. If None, the plot will be displayed.

    Returns:
    None
    """
    plt.rcParams.update({"font.size": FONTSIZE})
    fig, ax = plt.subplots(layout="constrained", dpi=DPI)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useOffset=False)
    lines = []
    for i, values in enumerate(plt_values):
        (line,) = ax.plot(time, values, label=labels[i])
        lines.append(line)
    ax.set_xlabel(xlabel)
    ax.grid()
    ax.legend()
    YLABEL = ", ".join(labels)
    ax.set_ylabel(YLABEL)
    if XLIM:
        ax.set_xlim(XLIM)
    if YLIM and not YZOOM:
        ax.set_ylim(YLIM)
    if YZOOM and not YLIM:
        if len(plt_values) > 1:
            print(f"YZOOM={YZOOM} is only allowed if one curve is plotted.")
            sys.exit(1)
        ax.set_ylim(plt_values[0][-1] - YZOOM / 2.0, plt_values[0][-1] + YZOOM / 2.0)
    if YLIM and YZOOM:
        print("YZOOM and YLIM not compatible. Choose one at set the other None.")
        sys.exit(1)

    if filename:
        if filename.endswith(".eps"):
            fig.savefig(filename, format="eps", transparent=TRANSPARENT)
        elif filename.endswith(".png"):
            fig.savefig(filename, format="png", transparent=TRANSPARENT)
        else:
            with open(filename, "w") as file:
                for i in range(len(time)):
                    file.write(f"{time[i]} {plt_values[i,:]}\n")
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
        "-s",
        "--step",
        action="store_true",
        help="use time steps rather than time in picoseconds",
    )
    parser.add_argument(
        "-u",
        "--unit",
        choices=["h", "ry", "ev", "kj/mol", "kcal/mol"],
        help="select energy unit",
    )
    parser.add_argument(
        "-e", "--etot", action="store_true", help="plot static total energy E_tot"
    )
    parser.add_argument(
        "-c",
        "--econ",
        action="store_true",
        help="plot conserved energy E_c in Hartree (or unit)",
    )
    parser.add_argument(
        "-f",
        "--fict",
        action="store_true",
        help="plot fictitious kinetic wave function energy in Hartree (or unit)",
    )
    parser.add_argument(
        "-t", "--temp", action="store_true", help="plot temperature in Kelvin"
    )
    parser.add_argument(
        "-a",
        choices=["p", "r"],
        action="append",
        help=(
            "plot friction\n" "type=p: for wave functions a_p\n" "type=r: for atoms a_r"
        ),
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
        "rootname", type=str, help="paw protocol file name without the .prot ending"
    )
    return parser.parse_args()


def get_unit(unit):
    """
    Converts a given unit to its corresponding tag and conversion factor.

    Parameters:
        unit (str): The unit to be converted.

    Returns:
        tuple: A tuple containing the unit tag and conversion factor.

    Example:
        >>> get_unit('ev')
        ('eV', 0.036749325)
    """
    unittag = "H"
    unitconv = 1.0
    if unit:
        unittag = unit.lower()
        if unittag == "h":
            unitconv = 1.0
            unittag = "H"
        elif unittag == "ry":
            unitconv = 0.5
            unittag = "Ry"
        elif unittag == "ev":
            unitconv = 1.0 / 27.211
            unittag = "eV"
        elif unittag == "kj/mol":
            unitconv = 1.0 / 2625.500223430069
            unittag = "kJ/mol"
        elif unittag == "kcal/mol":
            unitconv = 1.0 / 627.5096463920391
            unittag = "kcal/mol"
    return unittag, unitconv


def get_values_labels(values, args):
    """
    Generate values and labels for plotting based on the given arguments.

    Parameters:
        values (numpy.ndarray): The input values.
        args (argparse.Namespace): The command line arguments.

    Returns:
        tuple: A tuple containing the values and labels for plotting.
    """
    unittag, unitconv = get_unit(args.unit)
    plt_values = []
    plt_labels = []

    if args.etot:
        plt_values.append(values[:, 4] * unitconv)
        plt_labels.append(f"$E_{{tot}}$ [{unittag}]")
    if args.econ:
        plt_values.append(values[:, 5] * unitconv)
        plt_labels.append(f"$E_c$ [{unittag}]")
    if args.fict:
        plt_values.append(values[:, 3] * unitconv)
        plt_labels.append(rf"$E_f$ [{unittag}]")
    if args.temp:
        plt_values.append(values[:, 2])
        plt_labels.append(r"$T$ [K]")
    if args.a:
        for a in args.a:
            if a == "p":
                plt_values.append(values[:, 6])
                plt_labels.append(r"$a_p$")
            elif a == "r":
                plt_values.append(values[:, 7])
                plt_labels.append(r"$a_r$")
    return plt_values, plt_labels


def read_values_to_array(filepath, prefix="!>"):
    """
    Reads values from a file and stores them in a NumPy array.
    Args:
        filepath (str): The path to the file.
        prefix (str, optional): The prefix used to identify lines containing values. Defaults to '!>'.
    Returns:
        numpy.ndarray: A NumPy array containing the values read from the file.
    Raises:
        FileNotFoundError: If the file specified by `filepath` is not found.
        IOError: If there is an error reading the file.
        ValueError: If there is an error processing a line in the file.
    """
    values = []
    try:
        with open(filepath, "r") as file:
            for line in file:
                if line.startswith(prefix):
                    # Remove the prefix and any leading/trailing whitespace
                    line_values = line[len(prefix) :].strip()
                    # Split the line into values and convert them to floats
                    row = [float(value) for value in line_values.split()]
                    values.append(row)
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        return None
    except IOError as e:
        print(f"Error reading file: {e}")
        return None
    except ValueError as e:
        print(f"Error processing line: {e}")
        return None

    # Convert the list of rows into a NumPy array
    return np.array(values)


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
    values = read_values_to_array(args.rootname)
    if values is None:
        sys.exit(1)

    if args.step:
        time = values[:, 0]
        XLABEL = "time steps"
    else:
        time = values[:, 1]
        XLABEL = r"$t$ [ps]"

    plt_values, plt_labels = get_values_labels(values, args)
    plot_data(time, plt_values, plt_labels, XLABEL, args.out)


if __name__ == "__main__":
    main()
