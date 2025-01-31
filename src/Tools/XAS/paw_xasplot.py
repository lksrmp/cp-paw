import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

class XASData:
    def __init__(self, nkpt, nspin):
        self.nkpt = nkpt
        self.nspin = nspin
        self.data = {}
        self.nocc = {}

    def add_matrix(self, k_point, spin, matrix):
        if not (0 <= k_point < self.nkpt):
            raise ValueError(f'Invalid k_point: {k_point}. Must be between 0 and {self.nkpt - 1}.')
        if not (0 <= spin < self.nspin):
            raise ValueError(f'Invalid spin: {spin}. Must be between 0 and {self.nspin - 1}.')
        if k_point not in self.data:
            self.data[k_point] = {}
        self.data[k_point][spin] = matrix

    def add_nocc(self, k_point, spin, nocc):
        if not (0 <= k_point < self.nkpt):
            raise ValueError(f'Invalid k_point: {k_point}. Must be between 0 and {self.nkpt - 1}.')
        if not (0 <= spin < self.nspin):
            raise ValueError(f'Invalid spin: {spin}. Must be between 0 and {self.nspin - 1}.')
        if k_point not in self.data:
            self.nocc[k_point] = {}
        self.nocc[k_point][spin] = nocc

    def get_matrix(self, k_point, spin):
        if not (0 <= k_point < self.nkpt):
            raise ValueError(f'Invalid k_point: {k_point}. Must be between 0 and {self.nkpt - 1}.')
        if not (0 <= spin < self.nspin):
            raise ValueError(f'Invalid spin: {spin}. Must be between 0 and {self.nspin - 1}.')
        return self.data.get(k_point, {}).get(spin, None)
    
    def get_nocc(self, k_point, spin):
        if not (0 <= k_point < self.nkpt):
            raise ValueError(f'Invalid k_point: {k_point}. Must be between 0 and {self.nkpt - 1}.')
        if not (0 <= spin < self.nspin):
            raise ValueError(f'Invalid spin: {spin}. Must be between 0 and {self.nspin - 1}.')
        return self.nocc.get(k_point, {}).get(spin, None)

def main():
    parser = argparse.ArgumentParser(description='Visualize overlap matrix')
    parser.add_argument('filename', nargs='?', help='the file to process')
    parser.add_argument('--output-dir', help='directory to save the plots')
    parser.add_argument('--row-interval', type=str, help='row interval to select submatrix (e.g., "0:10")')
    parser.add_argument('--col-interval', type=str, help='column interval to select submatrix (e.g., "0:10")')
    args = parser.parse_args()

    if args.filename:
        XASData = read_data(args.filename)
    else:
        parser.print_help()
        sys.exit(1)

    row_interval = parse_interval(args.row_interval)
    col_interval = parse_interval(args.col_interval)

    for ikpt in range(1): #range(XASData.nkpt):
        for ispin in range(XASData.nspin):
            print(80*'#')
            print(f'## Processing k-point {ikpt+1}, spin {ispin+1}')
            print(80*'#')
            # plot_squared_eigenvalues_histogram(XASData, ikpt, ispin, args.output_dir)

            determinant = calculate_determinant(XASData, ikpt, ispin, occ=True)
            print(f'Determinant for k-point {ikpt+1}, spin {ispin+1}: {determinant}')

            # plot_row_vector_lengths(XASData, ikpt, ispin, args.output_dir)
            
            # plot_column_vector_lengths(XASData, ikpt, ispin, args.output_dir)

            plot_matrix(XASData, ikpt, ispin, args.output_dir, row_interval, col_interval)

            # plot_max_abs_values(XASData, ikpt, ispin, args.output_dir)

            # plot_max_values(XASData, ikpt, ispin, args.output_dir)

            # plot_linearly_dependent_columns(XASData, ikpt, ispin, args.output_dir)

            # plot_linearly_dependent_rows(XASData, ikpt, ispin, args.output_dir)

            cond = condition_number(XASData, ikpt, ispin, occ=False)
            print(f'Condition number for k-point {ikpt+1}, spin {ispin+1}: {cond}')

            # gram_matrix(XASData, ikpt, ispin)
            problem_rows = plot_singular_values(XASData, ikpt, ispin, args.output_dir, occ=True)
            print(f'Problematic rows: {problem_rows}')
            # if problem_rows:
            #     for row in problem_rows:
            #         exchange_determinant(XASData, ikpt, ispin, row, args.output_dir)
            


    if(args.output_dir):
        print(f'Plots saved to {args.output_dir}')
    else:
        plt.show()


def parse_interval(interval_str):
    if interval_str:
        try:
            start, end = map(int, interval_str.split(':'))
            return slice(start, end)
        except ValueError:
            print(f'Invalid interval format: {interval_str}. Expected format "start:end".')
            sys.exit(1)
    return slice(None)


def read_data(filename):
    # Placeholder for file processing logic
    print(f'Processing file: {filename}')
    with open(filename, 'rb') as f:
        # Read NKPT and NSPIN
        f.read(4)
        NKPT = np.fromfile(f, dtype=np.int32, count=1)[0]
        f.read(4)
        f.read(4)
        NSPIN = np.fromfile(f, dtype=np.int32, count=1)[0]
        f.read(4)
        print(NKPT, NSPIN)
        XAS = XASData(NKPT, NSPIN)
        for ikpt in range(NKPT):
            for ispin in range(NSPIN):
                f.read(4)
                nb2, nb1, nocc = np.fromfile(f, dtype=np.int32, count=3)
                f.read(4)
                print(ikpt, ispin, nb1, nb2, nocc)
                XAS.add_nocc(ikpt, ispin, nocc)
                f.read(4)
                data = np.fromfile(f, dtype=np.float64, count=2*nb1*nb2)
                f.read(4)
                matrix = data[0::2] + 1j * data[1::2]
                matrix = matrix.reshape((nb2, nb1), order='F')  # Use column-major order
                XAS.add_matrix(ikpt, ispin, matrix)
    return XAS

def plot_matrix(XASData, k_point, spin, output_dir=None, row_interval=slice(None), col_interval=slice(None)):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None:
        submatrix = matrix[row_interval, col_interval]
        plt.figure(dpi=150)
        plt.imshow(np.abs(submatrix), cmap='Blues', interpolation='none', vmin=0, extent=[col_interval.start or 0, col_interval.stop or matrix.shape[1], row_interval.stop or matrix.shape[0], row_interval.start or 0])
        plt.colorbar(label='abs. val.')
        # plt.imshow(np.abs(submatrix), cmap='Blues', interpolation='none', norm=LogNorm(vmin=1e-4, vmax=1.0), extent=[col_interval.start or 0, col_interval.stop or matrix.shape[1], row_interval.stop or matrix.shape[0], row_interval.start or 0])
        # plt.colorbar(label='abs. val. (log scale)')
        plt.title(f'Overlap Matrix for k-point {k_point+1}, spin {spin+1}')
        plt.xlabel('Column index')
        plt.ylabel('Row index')
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        # Draw a red border around the square sub matrix of size NOCC
        if row_interval.start:
            lenrow = nocc - row_interval.start
        else:
            lenrow = nocc+0.5
        if col_interval.start:
            lencol = nocc - col_interval.start
        else:
            lencol = nocc+0.5
        rect = plt.Rectangle((col_interval.start or 0 - 0.5, row_interval.start or 0 - 0.5), lencol, lenrow, edgecolor='red', facecolor='none', linewidth=2)
        plt.gca().add_patch(rect)
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'overlap_matrix_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')

def calculate_determinant(XASData, k_point, spin, row_interval=slice(None), col_interval=slice(None), occ=True):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None and nocc is not None:
        if occ:
            submatrix = matrix[:nocc, :nocc]
        else:
            submatrix = matrix
        determinant = np.linalg.det(submatrix)
        # # Check for near-zero rows or columns
        # row_sums = np.sum(np.abs(submatrix), axis=1)
        # col_sums = np.sum(np.abs(submatrix), axis=0)
        # # Find the lowest 10 row sums
        # lowest_row_sums_indices = np.argsort(row_sums)[:10]
        # lowest_row_sums = row_sums[lowest_row_sums_indices]
        # print(f'Lowest 10 row sums: {lowest_row_sums} at indices {lowest_row_sums_indices}')
        # # Find the lowest 10 column sums
        # lowest_col_sums_indices = np.argsort(col_sums)[:10]
        # lowest_col_sums = col_sums[lowest_col_sums_indices]
        # print(f'Lowest 10 column sums: {lowest_col_sums} at indices {lowest_col_sums_indices}')
        return determinant
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
        return None
        
def plot_squared_eigenvalues_histogram(XASData, k_point, spin, output_dir=None):
    matrix = XASData.get_matrix(k_point, spin)
    if matrix is not None:
        matrix = np.dot(matrix, matrix.conj().T)
        eigenvalues = np.linalg.eigvals(matrix)
        plt.figure(dpi=150)
        plt.hist(eigenvalues, bins=50, range=(eigenvalues.min(), eigenvalues.max()), alpha=0.75, label='Eigenvalues')
        plt.xlabel('Squared Eigenvalue')
        plt.ylabel('Frequency')
        plt.yscale('log')
        plt.title(f'Squared Eigenvalues Histogram for k-point {k_point+1}, spin {spin+1}')
        plt.legend()
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'squared_eigenvalues_histogram_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')

def plot_max_abs_values(XASData, k_point, spin, output_dir=None):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None:
        submatrix = matrix[:nocc, :nocc]
        max_abs_values = np.max(np.abs(submatrix), axis=0)
        plt.figure(dpi=150)
        plt.plot(max_abs_values, marker='o', linestyle='-', color='b')
        plt.xlabel('Column index')
        plt.ylabel('Max absolute value')
        plt.title(f'Max Absolute Values for k-point {k_point+1}, spin {spin+1}')
        plt.grid(True)
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'max_abs_values_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
        
def plot_max_values(XASData, k_point, spin, output_dir=None):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None:
        submatrix = matrix #[:nocc, :nocc]
        lenrow, lencol = submatrix.shape
        max_values_col = np.max(np.abs(submatrix), axis=0)
        max_positions_col = np.argmax(np.abs(submatrix), axis=0)
        max_values_row = np.max(np.abs(submatrix), axis=1)
        max_positions_row = np.argmax(np.abs(submatrix), axis=1)
        plt.figure(dpi=150)
        plt.scatter(max_positions_col, np.arange(0, lenrow), c=max_values_col, cmap='viridis', marker='o')
        plt.colorbar(label='Max absolute value')
        plt.xlabel('Column index')
        plt.ylabel('Row index')
        plt.title(f'Max Absolute Values and Positions for k-point {k_point+1}, spin {spin+1}')
        plt.grid(True)
        plt.figure(dpi=150)
        plt.scatter(np.arange(0, lencol), max_positions_row, c=max_values_row, cmap='viridis', marker='o')
        plt.colorbar(label='Max absolute value')
        plt.xlabel('Column index')
        plt.ylabel('Row index')
        plt.title(f'Max Absolute Values and Positions for k-point {k_point+1}, spin {spin+1}')
        plt.grid(True)
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'max_values_positions_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')

def find_linearly_dependent_columns(matrix, tol=1e-3):
    _, s, vh = np.linalg.svd(matrix)
    rank = np.sum(s > tol)
    dependent_columns = np.where(np.abs(s) <= tol)[0]
    return dependent_columns

def find_linearly_dependent_rows(matrix, tol=1e-3):
    return find_linearly_dependent_columns(matrix.T)

def plot_linearly_dependent_columns(XASData, k_point, spin, output_dir=None):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    submatrix = matrix[:nocc, :nocc]
    if matrix is not None:
        dependent_columns = find_linearly_dependent_columns(submatrix)
        plt.figure(dpi=150)
        plt.imshow(np.abs(submatrix), cmap='Blues', interpolation='none')
        plt.colorbar(label='abs. val.')
        plt.title(f'Linearly Dependent Columns for k-point {k_point+1}, spin {spin+1}')
        for col in dependent_columns:
            plt.axvline(x=col, color='red', linestyle='--')
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'linearly_dependent_columns_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')

def calculate_row_vector_lengths(XASData, k_point, spin):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None and nocc is not None:
        lengths = []
        for row_index in range(matrix.shape[0]):
            row_vector = matrix[row_index, :nocc]
            length = np.linalg.norm(row_vector)
            lengths.append(length)
        return lengths
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
        return None
    
def calculate_column_vector_lengths(XASData, k_point, spin):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None and nocc is not None:
        lengths = []
        for col_index in range(matrix.shape[1]):
            col_vector = matrix[:nocc, col_index]
            length = np.linalg.norm(col_vector)
            lengths.append(length)
        lengths = np.array(lengths)
        mult=1.0
        for i in range(nocc):
            mult=mult*lengths[i]
        print('Multiplied: '+str(mult))
        return lengths
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
        return None
    
def plot_row_vector_lengths(XASData, k_point, spin, output_dir=None):
    lengths = calculate_row_vector_lengths(XASData, k_point, spin)
    if lengths is not None:
        nocc = XASData.get_nocc(k_point, spin)
        plt.figure(dpi=150)
        plt.plot(lengths, marker='o', linestyle='-', color='b')
        plt.axvline(x=nocc-0.5, color='red', linestyle='--')
        plt.xlabel('Row index')
        plt.ylabel('Vector length')
        plt.title(f'Row Vector Lengths for k-point {k_point+1}, spin {spin+1}')
        plt.grid(True)
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'row_vector_lengths_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')

def plot_column_vector_lengths(XASData, k_point, spin, output_dir=None):
    lengths = calculate_column_vector_lengths(XASData, k_point, spin)
    if lengths is not None:
        nocc = XASData.get_nocc(k_point, spin)
        plt.figure(dpi=150)
        plt.plot(lengths, marker='o', linestyle='-', color='b')
        plt.axvline(x=nocc-0.5, color='red', linestyle='--')
        plt.xlabel('Column index')
        plt.ylabel('Vector length')
        plt.title(f'Column Vector Lengths for k-point {k_point+1}, spin {spin+1}')
        plt.grid(True)
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'column_vector_lengths_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
    
def plot_linearly_dependent_rows(XASData, k_point, spin, output_dir=None):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    submatrix = matrix[:nocc, :nocc]
    if matrix is not None:
        dependent_rows = find_linearly_dependent_rows(submatrix)
        plt.figure(dpi=150)
        plt.imshow(np.abs(submatrix), cmap='Blues', interpolation='none')
        plt.colorbar(label='abs. val.')
        plt.title(f'Linearly Dependent Rows for k-point {k_point+1}, spin {spin+1}')
        for row in dependent_rows:
            plt.axhline(y=row, color='red', linestyle='--')
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'linearly_dependent_rows_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')

def condition_number(XASData, k_point, spin, occ=True):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None and nocc is not None:
        if occ:
            submatrix = matrix[:nocc, :nocc]
        else:
            submatrix = matrix
        cond = np.linalg.cond(submatrix)
        return cond
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
        return None
    
def plot_singular_values(XASData, k_point, spin, output_dir=None, svdtol=0.5, conttol=0.6, occ=True):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None and nocc is not None:
        if occ:
            submatrix = matrix[:nocc, :nocc]
        else:
            submatrix = matrix
        U, singular_values, Vh = np.linalg.svd(submatrix)
        small_singular_values = np.where(singular_values < svdtol)[0]
        problematic_rows_return = []
        if len(small_singular_values) > 0:
            print(f"\nSmall singular values detected (threshold {svdtol}) (careful, Python indexing):")
            for idx in small_singular_values:
                print(f"Singular value {idx} is small: {singular_values[idx]}")
                # print(f"Corresponding right singular vector (v_{idx + 1}):")
                # print(Vh[idx])
                row_contributions = np.abs(Vh[idx])
                problematic_rows = np.where(row_contributions > conttol)[0]  # Adjust this threshold if needed
                problematic_rows_return.append(problematic_rows)
                print("Most influential rows (likely responsible for small singular value):", problematic_rows)
                print("Row contributions (magnitudes):", row_contributions[problematic_rows])
        else:
            print("\nNo small singular values detected.")
            problematic_rows_return = None
        determinant = np.prod(singular_values)
        print(f'Determinant from singular values for k-point {k_point+1}, spin {spin+1}: {determinant}')
        plt.figure(dpi=150)
        plt.plot(singular_values, marker='o', linestyle='-', color='b')
        plt.xlabel('Index')
        plt.ylabel('Singular value')
        plt.title(f'Singular Values for k-point {k_point+1}, spin {spin+1}')
        plt.grid(True)
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'singular_values_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
    return problematic_rows_return
    
def gram_matrix(XASData, k_point, spin):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None and nocc is not None:
        submatrix = matrix[:nocc, :nocc]
        gram = np.dot(submatrix, submatrix.conj().T)
        eigv = np.linalg.eigvals(gram)
        eigv = np.abs(eigv)
        plt.figure(dpi=150)
        plt.hist(eigv, bins=50, range=(eigv.min(), eigv.max()), alpha=0.75, label='Eigenvalues')
        plt.xlabel('Eigenvalue')
        plt.ylabel('Frequency')
        plt.yscale('log')
        plt.title(f'Gram Matrix Eigenvalues Histogram for k-point {k_point+1}, spin {spin+1}')
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
        return None
    
def exchange_rows_and_calculate_determinant(XASData, k_point, spin, row1, row2):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None and nocc is not None:
        submatrix = matrix[:nocc, :nocc].copy()
        if row1 < nocc and row2 >= nocc:
            submatrix[row1, :] = matrix[row2, :nocc]
            determinant = np.linalg.det(submatrix)
            return determinant
        else:
            print(f'Invalid row indices: row1={row1}, row2={row2}. Ensure row1 is occupied and row2 is unoccupied.')
            return None
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')
        return None
    
def exchange_determinant(XASData, k_point, spin, row1, output_dir=None):
    nocc = XASData.get_nocc(k_point, spin)
    nmax = XASData.get_matrix(k_point, spin).shape[0]
    det = []
    for i in range(nocc, nmax):
        det.append(exchange_rows_and_calculate_determinant(XASData, k_point, spin, row1, i))
    det = np.array(det)
    plt.figure(dpi=150)
    plt.plot(np.arange(nocc, nocc + len(det)), det, marker='o', linestyle='-', color='b')
    plt.xlabel('Row index')
    plt.ylabel('Determinant')
    plt.title(f'Determinant after exchanging row {row1} for k-point {k_point+1}, spin {spin+1}')
    plt.grid(True)
    if output_dir:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(os.path.join(output_dir, f'determinant_exchange_k{k_point+1}_s{spin+1}.png'))



if __name__ == '__main__':
    main()