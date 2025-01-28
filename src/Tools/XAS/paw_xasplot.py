import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

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
    args = parser.parse_args()

    if args.filename:
        XASData = read_data(args.filename)
    else:
        parser.print_help()
        sys.exit(1)

    for ikpt in range(XASData.nkpt):
        for ispin in range(XASData.nspin):
            plot_matrix(XASData, ikpt, ispin, args.output_dir)
    if(args.output_dir):
        print(f'Plots saved to {args.output_dir}')
    else:
        plt.show()
    # np.set_printoptions(precision=4, suppress=True, linewidth=np.inf)
    # print("Matrix (first 6x6 block):\n", matrix[0:6, 0:6])
    

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

def plot_matrix(XASData, k_point, spin, output_dir=None):
    matrix = XASData.get_matrix(k_point, spin)
    nocc = XASData.get_nocc(k_point, spin)
    if matrix is not None:
        plt.figure(dpi=150)
        plt.imshow(np.abs(matrix), cmap='Blues', interpolation='none', vmin=0)
        plt.colorbar(label='abs. val.')
        plt.title(f'Overlap Matrix for k-point {k_point+1}, spin {spin+1}')
        plt.xlabel('Column index')
        plt.ylabel('Row index')
        # Draw a red border around the square sub matrix of size NOCC
        rect = plt.Rectangle((-0.5, -0.5), nocc, nocc, edgecolor='red', facecolor='none', linewidth=2)
        plt.gca().add_patch(rect)
        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, f'overlap_matrix_k{k_point+1}_s{spin+1}.png'))
    else:
        print(f'No data available for k_point {k_point+1}, spin {spin+1}')

if __name__ == '__main__':
    main()