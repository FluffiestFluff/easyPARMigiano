###################################################################################################################
#   Automated Force Fields for Metals     /$$$$$$$   /$$$$$$  /$$$$$$$  /$$      /$$                              # 
#                                        | $$__  $$ /$$__  $$| $$__  $$| $$$    /$$$                              #
#   /$$$$$$   /$$$$$$   /$$$$$$$ /$$   /$$| $$  \ $$| $$  \ $$| $$  \ $$| $$$$  /$$$$                             #
#  /$$__  $$ |____  $$ /$$_____/| $$  | $$| $$$$$$$/| $$$$$$$$| $$$$$$$/| $$ $$/$$ $$                             #
# | $$$$$$$$  /$$$$$$$|  $$$$$$ | $$  | $$| $$____/ | $$__  $$| $$__  $$| $$  $$$| $$                             #
# | $$_____/ /$$__  $$ \____  $$| $$  | $$| $$      | $$  | $$| $$  \ $$| $$\  $ | $$                             #
# |  $$$$$$$|  $$$$$$$ /$$$$$$$/|  $$$$$$$| $$      | $$  | $$| $$  | $$| $$ \/  | $$                             #
#  \_______/ \_______/|_______/  \____  $$|__/      |__/  |__/|__/  |__/|__/     |__/                             #
#                               /$$  | $$                                                                         #
#                              |  $$$$$$/              Ver. 4.15 - 17 October 2025                                #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

import numpy as np
from collections import defaultdict
import sys
import periodictable
from functools import lru_cache

# function for gaussian input
def parse_gau(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None

    hessian_lines = []
    coord_lines = []
    charge_lines = []
    start_parsing_hessian = start_parsing_coords = start_parsing_charges = False
    coord_line_counter = 0
    
    for line in lines:
        # Parse Hessian
        if "Force constants in Cartesian coordinates" in line:
            start_parsing_hessian = True
            continue
        elif start_parsing_hessian and "Final forces over" in line:
            start_parsing_hessian = False
            continue
        elif start_parsing_hessian:
            hessian_lines.append(line.strip())
        
        # Parse coordinates
        if "Input orientation:" in line or "Z-Matrix orientation" in line:
            start_parsing_coords = True
            coord_line_counter = 0
            coord_lines = []
            continue
        elif start_parsing_coords:
            coord_line_counter += 1
            if coord_line_counter > 4:
                if line.strip().startswith("-----"):
                    start_parsing_coords = False
                else:
                    coord_lines.append(line.strip())
        
        # Parse Mulliken charges
        if "Mulliken charges:" in line or "Mulliken atomic charges:" in line:
            start_parsing_charges = True
            charge_lines = []
            continue
        elif start_parsing_charges and ("Sum of Mulliken charges" in line or "Sum of Mulliken atomic charges" in line):
            start_parsing_charges = False
        elif start_parsing_charges:
            charge_lines.append(line.strip())
    
    coordinates, atomic_numbers = extract_coordinates(coord_lines)
    hessian_size = 3 * len(atomic_numbers)
    hessian_matrix = parse_blocked_hessian(hessian_lines, hessian_size)
    charges = extract_charges(charge_lines)
    
    return coordinates, hessian_matrix, atomic_numbers, charges

# function for checkpoint input
def parse_fchk(filename):
    def parse_array(lines, startline, endline):
        arr = []
        in_section = False
        for line in lines:
            if line.startswith(startline):
                in_section = True
                continue
            if in_section and line.startswith(endline):
                break
            if in_section:
                try:
                    arr.extend(map(float, line.split()))
                except ValueError:
                    break
        return np.array(arr, dtype=float)

    with open(filename, "r") as f:
        content = f.readlines()

    # Try different end markers for coordinates
    end_markers = ['Number of symbols in', 'Force Field', 'Atomic numbers']
    for end_marker in end_markers:
        crds = parse_array(content, 'Current cartesian coordinates', end_marker)
        if len(crds) > 0:
            break

    # Try different end markers for Hessian
    end_markers = ['Nonadiabatic coupling', 'Dipole Moment', 'Vibrational Anharmonic']
    for end_marker in end_markers:
        hess_lower = parse_array(content, 'Cartesian Force Constants', end_marker)
        if len(hess_lower) > 0:
            break

    atomic_numbers = parse_array(content, 'Atomic numbers', 'Nuclear charges')
    charges = parse_array(content, 'Mulliken Charges', 'Cartesian Gradient')
 
    if len(crds) == 0 or len(hess_lower) == 0 or len(atomic_numbers) == 0:
        raise ValueError("Failed to parse required data from the file.")

    natoms = len(crds) // 3
    dim = 3 * natoms
    hess = np.zeros((dim, dim))
    idx = 0
    for i in range(dim):
        for j in range(i+1):
            hess[i, j] = hess_lower[idx]
            hess[j, i] = hess_lower[idx]
            idx += 1
    
    return crds.reshape(-1, 3), hess, atomic_numbers.astype(int), charges

# Helper functions for parse_gau
def extract_coordinates(coord_lines):
    coordinates = []
    atomic_numbers = []
    
    for line in coord_lines:
        parts = line.split()
        if len(parts) >= 6:
            atomic_numbers.append(int(parts[1]))
            coordinates.append(list(map(float, parts[3:6])))
    
    return np.array(coordinates), np.array(atomic_numbers)

def extract_charges(charge_lines):
    charges = []
    for line in charge_lines:
        parts = line.split()
        if len(parts) >= 2:
            try:
                charges.append(float(parts[-1]))
            except ValueError:
                continue
    return np.array(charges)

def parse_blocked_hessian(hessian_lines, hessian_size):
    hessian_matrix = np.zeros((hessian_size, hessian_size))
    current_col_start = 0
    
    for line in hessian_lines:
        parts = line.split()
        
        if all(p.isdigit() for p in parts):
            current_col_start = int(parts[0]) - 1
            continue
        
        row = int(parts[0]) - 1
        values = [float(val.replace('D', 'E')) for val in parts[1:]]
        
        for i, value in enumerate(values):
            col = current_col_start + i
            if col < hessian_size and row < hessian_size:
                hessian_matrix[row, col] = value
                hessian_matrix[col, row] = value
    
    return hessian_matrix

BOHR_TO_ANGSTROM = 0.52917720859  # Java FChkLoader constant
HARTREE_TO_KCAL_MOL = 627.509474
JAVA_ANGLE_FACTOR = 1971.3792352578043  # Java value: pi * Hartree_to_kcal


# Calculate distance between two coordinates.
def calculate_distance(coord1, coord2, is_fchk=False):
    """
    Return distance in Angstrom.

    Gaussian .fchk stores current cartesian coordinates in bohr.
    Gaussian output orientation coordinates are read here in Angstrom.
    """
    distance = np.linalg.norm(coord1 - coord2)
    if is_fchk:
        return distance * BOHR_TO_ANGSTROM
    return distance


def extract_sub_hessian(hessian, i, j):
    """
    Directed 3x3 Cartesian Hessian coupling block, Java/XYZViewer style.

    Java Hessian.makeSubMatrices() constructs submatrices[i][j] directly from
    the full Hessian data block:

        data[i*3 : i*3+3, j*3 : j*3+3]

    There is no leading minus sign in the Java code. Therefore ep2 returns
    H_ij directly, not -H_ij.
    """
    rows = [3 * i, 3 * i + 1, 3 * i + 2]
    cols = [3 * j, 3 * j + 1, 3 * j + 2]
    return hessian[np.ix_(rows, cols)]


def _java_cdiv(xr, xi, yr, yi):
    """
    Complex scalar division used by JAMA/Java EigenvalueDecomposition.

    Returns (qr, qi) = (xr + i xi) / (yr + i yi).
    """
    if abs(yr) > abs(yi):
        r = yi / yr
        d = yr + r * yi
        qr = (xr + r * xi) / d
        qi = (xi - r * xr) / d
    else:
        r = yr / yi
        d = yi + r * yr
        qr = (r * xr + xi) / d
        qi = (r * xi - xr) / d
    return qr, qi


def java_jama_eig_real(matrix):
    """
    Java/JAMA-style eigenvalue decomposition for a real nonsymmetric 3x3 matrix.

    This ports the nonsymmetric path used by the Java XYZViewer code:
        orthes()
        hqr2()
        cdiv()

    It returns the same kind of objects used by Matrix3D.getRealEigenvalues()
    and Matrix3D.getEigenvectors():

        d : real eigenvalue parts
        e : imaginary eigenvalue parts
        V : real eigenvector matrix columns

    Notes
    -----
    This is intentionally not np.linalg.eig(). NumPy returns complex
    eigenvectors for complex conjugate pairs, whereas Java/JAMA stores a real
    eigenvector matrix V and separate real/imaginary eigenvalue arrays. The
    Hessian force routines use only d[i] and real normalized columns V[:, i].
    """
    H = np.array(matrix, dtype=float, copy=True)
    n = H.shape[0]
    if H.shape != (n, n):
        raise ValueError("matrix must be square")
    if n != 3:
        raise ValueError("java_jama_eig_real is intended for 3x3 Hessian blocks")

    V = np.eye(n, dtype=float)
    d = np.zeros(n, dtype=float)
    e = np.zeros(n, dtype=float)
    ort = np.zeros(n, dtype=float)

    low = 0
    high = n - 1
    eps = 2.0 ** -52.0

    # orthes(): nonsymmetric reduction to Hessenberg form
    for m in range(low + 1, high):
        scale = 0.0
        for i in range(m, high + 1):
            scale += abs(H[i, m - 1])

        if scale != 0.0:
            h = 0.0
            for i in range(high, m - 1, -1):
                ort[i] = H[i, m - 1] / scale
                h += ort[i] * ort[i]

            g = np.sqrt(h)
            if ort[m] > 0.0:
                g = -g

            h -= ort[m] * g
            ort[m] -= g

            for j in range(m, n):
                f = 0.0
                for i in range(high, m - 1, -1):
                    f += ort[i] * H[i, j]
                f /= h
                for i in range(m, high + 1):
                    H[i, j] -= f * ort[i]

            for i in range(0, high + 1):
                f = 0.0
                for j in range(high, m - 1, -1):
                    f += ort[j] * H[i, j]
                f /= h
                for j in range(m, high + 1):
                    H[i, j] -= f * ort[j]

            ort[m] *= scale
            H[m, m - 1] = scale * g

    # Accumulate transformations.
    for m in range(high - 1, low, -1):
        if H[m, m - 1] != 0.0:
            for i in range(m + 1, high + 1):
                ort[i] = H[i, m - 1]
            for j in range(m, high + 1):
                g = 0.0
                for i in range(m, high + 1):
                    g += ort[i] * V[i, j]
                g = g / ort[m] / H[m, m - 1]
                for i in range(m, high + 1):
                    V[i, j] += g * ort[i]

    # hqr2(): nonsymmetric reduction from Hessenberg to real Schur form
    nn = n
    nidx = nn - 1
    exshift = 0.0
    p = q = r = s = z = 0.0

    norm = 0.0
    for i in range(nn):
        if (i < low) or (i > high):
            d[i] = H[i, i]
            e[i] = 0.0
        for j in range(max(i - 1, 0), nn):
            norm += abs(H[i, j])

    iter_count = 0
    while nidx >= low:
        l = nidx
        while l > low:
            s = abs(H[l - 1, l - 1]) + abs(H[l, l])
            if s == 0.0:
                s = norm
            if abs(H[l, l - 1]) < eps * s:
                break
            l -= 1

        # One root found.
        if l == nidx:
            H[nidx, nidx] += exshift
            d[nidx] = H[nidx, nidx]
            e[nidx] = 0.0
            nidx -= 1
            iter_count = 0
            continue

        # Two roots found.
        if l == nidx - 1:
            w = H[nidx, nidx - 1] * H[nidx - 1, nidx]
            p = (H[nidx - 1, nidx - 1] - H[nidx, nidx]) / 2.0
            q = p * p + w
            z = np.sqrt(abs(q))
            H[nidx, nidx] += exshift
            H[nidx - 1, nidx - 1] += exshift
            x = H[nidx, nidx]

            if q >= 0.0:
                if p >= 0.0:
                    z = p + z
                else:
                    z = p - z

                d[nidx - 1] = x + z
                d[nidx] = d[nidx - 1]
                if z != 0.0:
                    d[nidx] = x - w / z

                e[nidx - 1] = 0.0
                e[nidx] = 0.0

                x = H[nidx, nidx - 1]
                s = abs(x) + abs(z)
                p = x / s
                q = z / s
                r = np.sqrt(p * p + q * q)
                p /= r
                q /= r

                for j in range(nidx - 1, nn):
                    z = H[nidx - 1, j]
                    H[nidx - 1, j] = q * z + p * H[nidx, j]
                    H[nidx, j] = q * H[nidx, j] - p * z

                for i in range(0, nidx + 1):
                    z = H[i, nidx - 1]
                    H[i, nidx - 1] = q * z + p * H[i, nidx]
                    H[i, nidx] = q * H[i, nidx] - p * z

                for i in range(low, high + 1):
                    z = V[i, nidx - 1]
                    V[i, nidx - 1] = q * z + p * V[i, nidx]
                    V[i, nidx] = q * V[i, nidx] - p * z
            else:
                d[nidx - 1] = x + p
                d[nidx] = x + p
                e[nidx - 1] = z
                e[nidx] = -z

            nidx -= 2
            iter_count = 0
            continue

        # No convergence yet.
        x = H[nidx, nidx]
        y = 0.0
        w = 0.0
        if l < nidx:
            y = H[nidx - 1, nidx - 1]
            w = H[nidx, nidx - 1] * H[nidx - 1, nidx]

        # Wilkinson's original ad hoc shift.
        if iter_count == 10:
            exshift += x
            for i in range(low, nidx + 1):
                H[i, i] -= x
            s = abs(H[nidx, nidx - 1]) + abs(H[nidx - 1, nidx - 2])
            x = y = 0.75 * s
            w = -0.4375 * s * s

        # MATLAB's new ad hoc shift.
        if iter_count == 30:
            s = (y - x) / 2.0
            s = s * s + w
            if s > 0.0:
                s = np.sqrt(s)
                if y < x:
                    s = -s
                s = x - w / ((y - x) / 2.0 + s)
                for i in range(low, nidx + 1):
                    H[i, i] -= s
                exshift += s
                x = y = w = 0.964

        iter_count += 1

        # Look for two consecutive small sub-diagonal elements.
        m = nidx - 2
        while m >= l:
            z = H[m, m]
            r = x - z
            s = y - z
            p = (r * s - w) / H[m + 1, m] + H[m, m + 1]
            q = H[m + 1, m + 1] - z - r - s
            r = H[m + 2, m + 1]
            s = abs(p) + abs(q) + abs(r)
            p /= s
            q /= s
            r /= s
            if m == l:
                break
            if abs(H[m, m - 1]) * (abs(q) + abs(r)) < (
                eps * abs(p) *
                (abs(H[m - 1, m - 1]) + abs(z) + abs(H[m + 1, m + 1]))
            ):
                break
            m -= 1

        for i in range(m + 2, nidx + 1):
            H[i, i - 2] = 0.0
            if i > m + 2:
                H[i, i - 3] = 0.0

        # Double QR step involving rows l:n and columns m:n.
        for k in range(m, nidx):
            notlast = (k != nidx - 1)
            if k != m:
                p = H[k, k - 1]
                q = H[k + 1, k - 1]
                r = H[k + 2, k - 1] if notlast else 0.0
                x = abs(p) + abs(q) + abs(r)
                if x != 0.0:
                    p /= x
                    q /= x
                    r /= x

            if x == 0.0:
                break

            s = np.sqrt(p * p + q * q + r * r)
            if p < 0.0:
                s = -s

            if s != 0.0:
                if k != m:
                    H[k, k - 1] = -s * x
                elif l != m:
                    H[k, k - 1] = -H[k, k - 1]

                p += s
                x = p / s
                y = q / s
                z = r / s
                q /= p
                r /= p

                for j in range(k, nn):
                    p = H[k, j] + q * H[k + 1, j]
                    if notlast:
                        p += r * H[k + 2, j]
                        H[k + 2, j] -= p * z
                    H[k, j] -= p * x
                    H[k + 1, j] -= p * y

                for i in range(0, min(nidx, k + 3) + 1):
                    p = x * H[i, k] + y * H[i, k + 1]
                    if notlast:
                        p += z * H[i, k + 2]
                        H[i, k + 2] -= p * r
                    H[i, k] -= p
                    H[i, k + 1] -= p * q

                for i in range(low, high + 1):
                    p = x * V[i, k] + y * V[i, k + 1]
                    if notlast:
                        p += z * V[i, k + 2]
                        V[i, k + 2] -= p * r
                    V[i, k] -= p
                    V[i, k + 1] -= p * q

    if norm == 0.0:
        return d, e, V

    # Backsubstitute to find vectors of upper triangular form.
    for nidx in range(nn - 1, -1, -1):
        p = d[nidx]
        q = e[nidx]

        if q == 0.0:
            l = nidx
            H[nidx, nidx] = 1.0

            for i in range(nidx - 1, -1, -1):
                w = H[i, i] - p
                r = 0.0
                for j in range(l, nidx + 1):
                    r += H[i, j] * H[j, nidx]

                if e[i] < 0.0:
                    z = w
                    s = r
                else:
                    l = i
                    if e[i] == 0.0:
                        if w != 0.0:
                            H[i, nidx] = -r / w
                        else:
                            H[i, nidx] = -r / eps * norm
                    else:
                        x = H[i, i + 1]
                        y = H[i + 1, i]
                        q2 = (d[i] - p) * (d[i] - p) + e[i] * e[i]
                        t = (x * s - z * r) / q2
                        H[i, nidx] = t
                        if abs(x) > abs(z):
                            H[i + 1, nidx] = (-r - w * t) / x
                        else:
                            H[i + 1, nidx] = (-s - y * t) / z

                    t = abs(H[i, nidx])
                    if eps * t * t > 1.0:
                        for j in range(i, nidx + 1):
                            H[j, nidx] /= t

        elif q < 0.0:
            l = nidx - 1

            if abs(H[nidx, nidx - 1]) > abs(H[nidx - 1, nidx]):
                H[nidx - 1, nidx - 1] = q / H[nidx, nidx - 1]
                H[nidx - 1, nidx] = -(H[nidx, nidx] - p) / H[nidx, nidx - 1]
            else:
                cr, ci = _java_cdiv(0.0, -H[nidx - 1, nidx], H[nidx - 1, nidx - 1] - p, q)
                H[nidx - 1, nidx - 1] = cr
                H[nidx - 1, nidx] = ci

            H[nidx, nidx - 1] = 0.0
            H[nidx, nidx] = 1.0

            for i in range(nidx - 2, -1, -1):
                ra = 0.0
                sa = 0.0
                for j in range(l, nidx + 1):
                    ra += H[i, j] * H[j, nidx - 1]
                    sa += H[i, j] * H[j, nidx]

                w = H[i, i] - p

                if e[i] < 0.0:
                    z = w
                    r = ra
                    s = sa
                else:
                    l = i
                    if e[i] == 0.0:
                        cr, ci = _java_cdiv(-ra, -sa, w, q)
                        H[i, nidx - 1] = cr
                        H[i, nidx] = ci
                    else:
                        x = H[i, i + 1]
                        y = H[i + 1, i]
                        vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q
                        vi = (d[i] - p) * 2.0 * q
                        if (vr == 0.0) and (vi == 0.0):
                            vr = eps * norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(z))

                        cr, ci = _java_cdiv(
                            x * r - z * ra + q * sa,
                            x * s - z * sa - q * ra,
                            vr,
                            vi,
                        )
                        H[i, nidx - 1] = cr
                        H[i, nidx] = ci

                        if abs(x) > abs(z) + abs(q):
                            H[i + 1, nidx - 1] = (-ra - w * H[i, nidx - 1] + q * H[i, nidx]) / x
                            H[i + 1, nidx] = (-sa - w * H[i, nidx] - q * H[i, nidx - 1]) / x
                        else:
                            cr, ci = _java_cdiv(
                                -r - y * H[i, nidx - 1],
                                -s - y * H[i, nidx],
                                z,
                                q,
                            )
                            H[i + 1, nidx - 1] = cr
                            H[i + 1, nidx] = ci

                    t = max(abs(H[i, nidx - 1]), abs(H[i, nidx]))
                    if eps * t * t > 1.0:
                        for j in range(i, nidx + 1):
                            H[j, nidx - 1] /= t
                            H[j, nidx] /= t

    # Back transformation to original basis.
    for j in range(nn - 1, low - 1, -1):
        for i in range(low, high + 1):
            z = 0.0
            for k in range(low, min(j, high) + 1):
                z += V[i, k] * H[k, j]
            V[i, j] = z

    return d, e, V


def java_matrix3d_real_eigenvectors(matrix, tol=1e-14):
    """
    Return Java Matrix3D.getRealEigenvalues() and getEigenvectors() equivalents.

    Matrix3D.getEigenvectors() normalizes each real column returned from the
    internal eigenVec matrix. It does not use complex NumPy eigenvectors.
    """
    d, e, V = java_jama_eig_real(matrix)
    vectors = []
    for i in range(3):
        vec = np.array([V[0, i], V[1, i], V[2, i]], dtype=float)
        nrm = np.linalg.norm(vec)
        if nrm < tol:
            vectors.append(vec)
        else:
            vectors.append(vec / nrm)
    return d, e, vectors


def _bond_force_one_direction(hessian, coordinates, idx1, idx2):
    """
    Java/XYZViewer one-direction Seminario bond force constant.

    Java Hessian.calcForce():
        Vector3D[] v = submatrices[i][j].getEigenvectors()
        double[] lambda = submatrices[i][j].getRealEigenvalues()
        kij += lambda[n] * abs(unit_bond.dot(v[n]))
        return 627.5095 * kij / 0.2800285180911043

    Here the eigenvectors/eigenvalues are generated by the Python port of the
    Java/JAMA nonsymmetric eigenvalue decomposition, not by np.linalg.eig().
    """
    sub_hessian = extract_sub_hessian(hessian, idx1, idx2)

    vec12 = coordinates[idx1] - coordinates[idx2]
    unit_vec12 = vec12 / np.linalg.norm(vec12)

    lambdas, _imag, vectors = java_matrix3d_real_eigenvectors(sub_hessian)

    kij = 0.0
    for lam, vec in zip(lambdas, vectors):
        kij += lam * abs(np.dot(unit_vec12, vec))

    # Java calcForce() returns the signed converted value.
    # The GUI-level final value is commonly abs(k12 + k21) / 2.
    return kij * HARTREE_TO_KCAL_MOL / (BOHR_TO_ANGSTROM ** 2)


def calculate_bond_force_constant(hessian, coordinates, idx1, idx2):
    """
    Java-like final bond constant.

    Java GUI reports the average of calcForce(i,j) and calcForce(j,i).
    """
    k12 = _bond_force_one_direction(hessian, coordinates, idx1, idx2)
    k21 = _bond_force_one_direction(hessian, coordinates, idx2, idx1)
    return abs(k12 + k21) / 2.0


def seminario_method(hessian, coordinates, atom_pairs, atomic_numbers):
    bonds = []
    for atom1, atom2 in atom_pairs:
        idx1, idx2 = atom1 - 1, atom2 - 1
        force_constant = calculate_bond_force_constant(hessian, coordinates, idx1, idx2)
        bonds.append((atom1, atom2, force_constant))
    return bonds


def read_distances(filename):
    atom_pairs = []
    with open(filename, "r") as f:
        for line in f:
            if line.strip():
                parts = line.split()
                atom1 = int(parts[0])
                atom2 = int(parts[1])
                atom_pairs.append((atom1, atom2))
    return atom_pairs


def read_angles(filename):
    angles = []
    with open(filename, "r") as f:
        for line in f:
            if line.strip():
                parts = line.split()
                atom1 = int(parts[0])
                atom2 = int(parts[1])
                atom3 = int(parts[2])
                angle_value = float(parts[3])
                angles.append((atom1, atom2, atom3, angle_value))
    return angles


def calculate_angle(coord1, coord2, coord3):
    """
    Return angle in degrees.

    The value is scale-independent, so bohr vs Angstrom does not matter.
    """
    v1 = coord1 - coord2
    v2 = coord3 - coord2

    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))

    return np.degrees(angle)


def _coords_for_java_geometry(coordinates, is_fchk=False):
    """
    Return coordinates in the same length convention as the Java viewer.

    Evidence from the Java MolData exporter: Turbomole output divides internal
    coordinates by 0.529177..., which means Java internal coordinates are Angstrom.
    easyPARM's parse_fchk returns Gaussian .fchk coordinates in bohr, so .fchk
    geometries must be converted to Angstrom before using Java-like angle lengths.
    """
    if is_fchk:
        return coordinates * BOHR_TO_ANGSTROM
    return coordinates


def _angle_force_one_direction(hessian, coordinates, idx1, idx2, idx3, is_fchk=False):
    """
    Java-like one-direction A-B-C angle force constant.

    Mirrors Java calcAngleForce(a,b,c,...):
        Uab = A - B
        Ucb = C - B
        Un  = Ucb x Uab
        Upa = Un x Uab
        Upc = Ucb x Un

        k_AB = sum(lambda_i * abs(dot(Upa, eigvec_i))) * R_AB^2
        k_CB = sum(lambda_i * abs(dot(Upc, eigvec_i))) * R_CB^2

        ktheta_inv = 1/k_AB + 1/k_CB
        ktheta     = JAVA_ANGLE_FACTOR / ktheta_inv

    This intentionally reproduces the Java convention, including the
    JAVA_ANGLE_FACTOR. It is Java-compatible, not a re-derivation of Amber units.
    """
    geom = _coords_for_java_geometry(coordinates, is_fchk=is_fchk)

    Uab = geom[idx1] - geom[idx2]
    Ucb = geom[idx3] - geom[idx2]

    Un = np.cross(Ucb, Uab)
    Un = Un / np.linalg.norm(Un)

    Upa = np.cross(Un, Uab)
    Upc = np.cross(Ucb, Un)

    Upa = Upa / np.linalg.norm(Upa)
    Upc = Upc / np.linalg.norm(Upc)

    Rab = np.linalg.norm(Uab)
    Rcb = np.linalg.norm(Ucb)

    # Java calcAngleForce(a,b,c) uses submatrices[b][a] and submatrices[b][c].
    sub_hessian_BA = extract_sub_hessian(hessian, idx2, idx1)
    sub_hessian_BC = extract_sub_hessian(hessian, idx2, idx3)

    lambdas_BA, _imag_BA, vectors_BA = java_matrix3d_real_eigenvectors(sub_hessian_BA)
    lambdas_BC, _imag_BC, vectors_BC = java_matrix3d_real_eigenvectors(sub_hessian_BC)

    k_ab = 0.0
    for lam, vec in zip(lambdas_BA, vectors_BA):
        k_ab += lam * abs(np.dot(Upa, vec))
    k_ab *= Rab * Rab

    k_cb = 0.0
    for lam, vec in zip(lambdas_BC, vectors_BC):
        k_cb += lam * abs(np.dot(Upc, vec))
    k_cb *= Rcb * Rcb

    ktheta_inv = (1.0 / k_ab) + (1.0 / k_cb)
    return JAVA_ANGLE_FACTOR / ktheta_inv


def calculate_angle_force_constant(hessian, coordinates, idx1, idx2, idx3, is_fchk=False):
    """
    Java-like final angle force constant.

    Java GUI reports the average of A-B-C and C-B-A calculations.
    """
    k_forward = _angle_force_one_direction(
        hessian, coordinates, idx1, idx2, idx3, is_fchk=is_fchk
    )
    k_reverse = _angle_force_one_direction(
        hessian, coordinates, idx3, idx2, idx1, is_fchk=is_fchk
    )
    return abs(k_forward + k_reverse) / 2.0


def read_dihedrals(filename):
    dihedrals = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.split()
                if len(parts) == 5:
                    atom1, atom2, atom3, atom4, dihedral_value = map(float, parts)
                    dihedrals.append((int(atom1), int(atom2), int(atom3), int(atom4), dihedral_value))
    return dihedrals

# Enhanced molecular connectivity analysis
def analyze_molecular_connectivity(coordinates, atomic_numbers, atom_names):
    
    bohr_to_angstrom = 0.529177
    natom = len(atomic_numbers)
    
    # Vectorized distance calculation using broadcasting
    coords_expanded_i = coordinates[:, np.newaxis, :]
    coords_expanded_j = coordinates[np.newaxis, :, :]
    distance_matrix = np.linalg.norm(coords_expanded_i - coords_expanded_j, axis=2) * bohr_to_angstrom
    np.fill_diagonal(distance_matrix, 0)  # Ensure diagonal is 0
    
    # Define typical bond lengths for connectivity (in Angstroms)
    bond_cutoffs = {
        ('H', 'C'): 1.2, ('H', 'N'): 1.2, ('H', 'O'): 1.2,
        ('C', 'C'): 1.8, ('C', 'N'): 1.8, ('C', 'O'): 1.8,
        ('N', 'N'): 1.8, ('N', 'O'): 1.8, ('O', 'O'): 1.8,
        ('H', 'S'): 1.5, ('C', 'S'): 2.0, ('S', 'S'): 2.3,
        ('H', 'P'): 1.5, ('C', 'P'): 2.0, ('P', 'P'): 2.3,
        # Metal bonds (more generous cutoffs)
        ('Metal', 'H'): 2.0, ('Metal', 'C'): 2.5, ('Metal', 'N'): 2.5,
        ('Metal', 'O'): 2.5, ('Metal', 'S'): 2.8, ('Metal', 'P'): 2.8,
        ('Metal', 'Metal'): 3.0
    }
    
    # Metal atoms set
    metal_atoms = {
        'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
        'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
        'LA', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
        'AC', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN'
    }
    
    def get_cutoff(elem1, elem2):
        elem1_is_metal = periodictable.elements[elem1].symbol.upper() in metal_atoms
        elem2_is_metal = periodictable.elements[elem2].symbol.upper() in metal_atoms
        
        if elem1_is_metal or elem2_is_metal:
            if elem1_is_metal and elem2_is_metal:
                return bond_cutoffs.get(('Metal', 'Metal'), 3.0)
            else:
                non_metal = elem2 if elem1_is_metal else elem1
                non_metal_symbol = periodictable.elements[non_metal].symbol.upper()
                return bond_cutoffs.get(('Metal', non_metal_symbol), 
                                      bond_cutoffs.get((non_metal_symbol, 'Metal'), 2.5))
        else:
            elem1_symbol = periodictable.elements[elem1].symbol.upper()
            elem2_symbol = periodictable.elements[elem2].symbol.upper()
            pair1 = (elem1_symbol, elem2_symbol)
            pair2 = (elem2_symbol, elem1_symbol)
            return bond_cutoffs.get(pair1, bond_cutoffs.get(pair2, 1.8))
    
    # Pre-compute cutoff matrix
    cutoff_matrix = np.zeros((natom, natom))
    for i in range(natom):
        for j in range(i+1, natom):
            cutoff = get_cutoff(atomic_numbers[i], atomic_numbers[j])
            cutoff_matrix[i, j] = cutoff
            cutoff_matrix[j, i] = cutoff
    
    # Vectorized connectivity determination
    connectivity = distance_matrix <= cutoff_matrix
    np.fill_diagonal(connectivity, False)
    
    # Create neighbor lists
    neighbors = {i: np.where(connectivity[i])[0].tolist() for i in range(natom)}
    
    return connectivity, neighbors, distance_matrix

#Create a detailed chemical fingerprint for an atom based on its environment
def get_chemical_environment_fingerprint(atomic_numbers, atom_idx, neighbors, max_depth=3, _cache=None):
    
    # Use memoization to avoid redundant calculations
    if _cache is None:
        _cache = {}
    
    cache_key = (atom_idx, max_depth)
    if cache_key in _cache:
        return _cache[cache_key]
    
    def explore_environment(current_atom, depth, visited, path):
        if depth > max_depth:
            return []
        
        visited.add(current_atom)
        environment_info = []
        
        # Get immediate neighbors
        current_neighbors = []
        for neighbor_idx in neighbors[current_atom]:
            if neighbor_idx not in visited:
                neighbor_element = periodictable.elements[atomic_numbers[neighbor_idx]].symbol.upper()
                current_neighbors.append((neighbor_element, neighbor_idx))
        
        # Sort neighbors for consistent fingerprinting
        current_neighbors.sort(key=lambda x: (x[0], x[1]))
        
        for neighbor_element, neighbor_idx in current_neighbors:
            # Create path signature
            new_path = path + [neighbor_element]
            environment_info.append(tuple(new_path))
            
            # Recursively explore deeper
            if depth < max_depth:
                deeper_env = explore_environment(neighbor_idx, depth + 1, visited.copy(), new_path)
                environment_info.extend(deeper_env)
        
        return environment_info
    
    # Start exploration from the atom
    fingerprint = tuple(sorted(explore_environment(atom_idx, 0, set(), [periodictable.elements[atomic_numbers[atom_idx]].symbol.upper()])))
    _cache[cache_key] = fingerprint
    return fingerprint

#Identify aromatic rings in the molecule using simple ring detection
def identify_aromatic_rings(atomic_numbers, connectivity, neighbors):
    
    natom = len(atomic_numbers)
    aromatic_atoms = set()
    
    # Use iterative approach with visited tracking to avoid redundant searches
    def find_rings_from_atom(start_atom, max_ring_size=8):
        stack = [(start_atom, [start_atom], set())]
        rings_found = []
        
        while stack:
            current, path, visited_edges = stack.pop()
            
            if len(path) > max_ring_size:
                continue
            
            for neighbor in neighbors[current]:
                edge = tuple(sorted([current, neighbor]))
                
                if neighbor == start_atom and len(path) >= 3:
                    # Found a ring
                    rings_found.append(path[:])
                elif neighbor not in path and edge not in visited_edges:
                    new_visited = visited_edges.copy()
                    new_visited.add(edge)
                    stack.append((neighbor, path + [neighbor], new_visited))
        
        return rings_found
    
    # Find all rings
    all_rings = []
    processed = set()
    
    for atom in range(natom):
        if atom not in processed:
            rings = find_rings_from_atom(atom)
            for ring in rings:
                if len(ring) >= 5 and len(ring) <= 8:  # Reasonable ring sizes
                    ring_tuple = tuple(sorted(ring))
                    if ring_tuple not in [tuple(sorted(r)) for r in all_rings]:
                        all_rings.append(ring)
                        processed.update(ring)
    
    # Check for aromaticity (simplified heuristic)
    for ring in all_rings:
        if len(ring) in [5, 6]:  # Common aromatic ring sizes
            # Check if ring is mostly carbons with some sp2 character
            carbon_count = sum(1 for atom in ring if periodictable.elements[atomic_numbers[atom]].symbol.upper() == 'C')
            nitrogen_count = sum(1 for atom in ring if periodictable.elements[atomic_numbers[atom]].symbol.upper() == 'N')
            
            # Heuristic: if mostly C and N, likely aromatic
            if (carbon_count + nitrogen_count) >= len(ring) * 0.8:
                aromatic_atoms.update(ring)
    
    return aromatic_atoms

#Estimate carbon hybridization based on bonding pattern.
def analyze_carbon_hybridization(atomic_numbers, carbon_idx, neighbors):
    carbon_neighbors = neighbors[carbon_idx]
    num_neighbors = len(carbon_neighbors)
    
    # Count different atom types
    neighbor_elements = [periodictable.elements[atomic_numbers[n]].symbol.upper() for n in carbon_neighbors]
    hydrogen_count = neighbor_elements.count('H')
    heavy_atom_count = num_neighbors - hydrogen_count
    
    if num_neighbors == 4:
        return 'sp3'
    elif num_neighbors == 3:
        return 'sp2'
    elif num_neighbors == 2:
        return 'sp'
    else:
        return 'unknown'

#Enhanced hydrogen environment analysis 
def analyze_enhanced_hydrogen_environments(atomic_numbers, neighbors, connectivity, aromatic_atoms):
    natom = len(atomic_numbers)
    
    # Categorize hydrogens with detailed environment analysis
    hydrogen_environments = {}
    
    # Pre-filter hydrogen atoms
    hydrogen_indices = [i for i in range(natom) if periodictable.elements[atomic_numbers[i]].symbol.upper() == 'H']
    
    for i in hydrogen_indices:
        # Get the heavy atom this H is bonded to
        heavy_neighbors = []
        for neighbor_idx in neighbors[i]:
            neighbor_element = periodictable.elements[atomic_numbers[neighbor_idx]].symbol.upper()
            if neighbor_element != 'H':
                heavy_neighbors.append((neighbor_idx, neighbor_element))
        
        if len(heavy_neighbors) != 1:
            continue  # Skip unusual cases
            
        heavy_atom_idx, heavy_element = heavy_neighbors[0]
        
        # Analyze specific environments
        if heavy_element == 'C':
            carbon_hybridization = analyze_carbon_hybridization(atomic_numbers, heavy_atom_idx, neighbors)
            is_aromatic_carbon = heavy_atom_idx in aromatic_atoms
            
            # Get carbon's neighbors (excluding this hydrogen)
            carbon_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            carbon_neighbor_elements = tuple(sorted([periodictable.elements[atomic_numbers[n]].symbol.upper() for n in carbon_neighbors]))
            
            # Count different types of neighbors
            h_on_carbon = carbon_neighbor_elements.count('H')
            c_on_carbon = carbon_neighbor_elements.count('C')
            
            # Create detailed environment key
            if is_aromatic_carbon:
                env_key = f"aromatic_H_{carbon_hybridization}_{carbon_neighbor_elements}"
            else:
                # Aliphatic carbon - be more specific about environment
                if h_on_carbon == 2 and c_on_carbon == 1:
                    # Check what the carbon neighbor is connected to
                    carbon_neighbor_idx = [n for n in carbon_neighbors if periodictable.elements[atomic_numbers[n]].symbol.upper() == 'C'][0]
                    carbon_neighbor_neighbors = tuple(sorted([periodictable.elements[atomic_numbers[n]].symbol.upper() for n in neighbors[carbon_neighbor_idx]]))
                    env_key = f"methyl_H_connected_to_C_with_{carbon_neighbor_neighbors}"
                elif h_on_carbon == 1 and c_on_carbon == 2:
                    # Methylene - check what carbons are connected to
                    carbon_neighbor_envs = []
                    for cn_idx in [n for n in carbon_neighbors if periodictable.elements[atomic_numbers[n]].symbol.upper() == 'C']:
                        cn_neighbors = tuple(sorted([periodictable.elements[atomic_numbers[n]].symbol.upper() for n in neighbors[cn_idx]]))
                        carbon_neighbor_envs.append(cn_neighbors)
                    env_key = f"methylene_H_{tuple(sorted(carbon_neighbor_envs))}"
                elif h_on_carbon == 0:
                    # Tertiary or quaternary carbon
                    env_key = f"tertiary_H_{carbon_neighbor_elements}"
                else:
                    # General case
                    env_key = f"aliphatic_H_{carbon_hybridization}_{carbon_neighbor_elements}"
                    
        elif heavy_element == 'N':
            # Nitrogen-bonded hydrogens
            nitrogen_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            nitrogen_neighbor_elements = tuple(sorted([periodictable.elements[atomic_numbers[n]].symbol.upper() for n in nitrogen_neighbors]))
            is_aromatic_nitrogen = heavy_atom_idx in aromatic_atoms
            
            if is_aromatic_nitrogen:
                env_key = f"aromatic_NH_{nitrogen_neighbor_elements}"
            else:
                env_key = f"amine_H_{nitrogen_neighbor_elements}"
                
        elif heavy_element == 'O':
            # Oxygen-bonded hydrogens
            oxygen_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            oxygen_neighbor_elements = tuple(sorted([periodictable.elements[atomic_numbers[n]].symbol.upper() for n in oxygen_neighbors]))
            env_key = f"hydroxyl_H_{oxygen_neighbor_elements}"
            
        elif heavy_element == 'S':
            # Sulfur-bonded hydrogens
            sulfur_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            sulfur_neighbor_elements = tuple(sorted([periodictable.elements[atomic_numbers[n]].symbol.upper() for n in sulfur_neighbors]))
            env_key = f"thiol_H_{sulfur_neighbor_elements}"
            
        else:
            # Other heavy atoms
            env_key = f"{heavy_element}_H"
        
        # Group hydrogens with identical environments
        if env_key not in hydrogen_environments:
            hydrogen_environments[env_key] = []
        hydrogen_environments[env_key].append(i)
    
    # Create constraint groups for hydrogens with identical environments
    hydrogen_constraint_groups = []
    for env_key, h_indices in hydrogen_environments.items():
        if len(h_indices) > 1:
            hydrogen_constraint_groups.append({
                'type': 'equivalent_hydrogens',
                'indices': h_indices,
                'description': f"Equivalent hydrogens: {env_key}",
                'environment': env_key
            })
    
    return hydrogen_constraint_groups

#Enhanced equivalent atom identification using chemical environment fingerprints
def identify_equivalent_atoms_enhanced(coordinates, atomic_numbers, atom_names):
    
    connectivity, neighbors, distance_matrix = analyze_molecular_connectivity(coordinates, atomic_numbers, atom_names)
    
    natom = len(atomic_numbers)
    equivalent_groups = []
    processed = set()
    
    # Pre-identify aromatic atoms once
    aromatic_atoms = identify_aromatic_rings(atomic_numbers, connectivity, neighbors)
    
    # Group atoms by element first
    elements_dict = defaultdict(list)
    for i in range(natom):
        element = periodictable.elements[atomic_numbers[i]].symbol.upper()
        elements_dict[element].append(i)
    
    # Use shared cache for fingerprint calculations
    fingerprint_cache = {}
    
    # For each element, find truly equivalent atoms
    for element, atom_indices in elements_dict.items():
        if len(atom_indices) < 2:
            continue
            
        # Create environment fingerprints for all atoms of this element
        fingerprints = defaultdict(list)
        for atom_idx in atom_indices:
            if atom_idx in processed:
                continue
                
            fingerprint = get_chemical_environment_fingerprint(atomic_numbers, atom_idx, neighbors, max_depth=3, _cache=fingerprint_cache)
            fingerprints[fingerprint].append(atom_idx)
        
        # Create equivalent groups for atoms with identical fingerprints
        for fingerprint, equivalent_atoms in fingerprints.items():
            if len(equivalent_atoms) > 1:
                equivalent_groups.append(equivalent_atoms)
                processed.update(equivalent_atoms)
                
    # Also analyze hydrogen environments specifically - pass pre-computed aromatic_atoms
    hydrogen_groups = analyze_enhanced_hydrogen_environments(atomic_numbers, neighbors, connectivity, aromatic_atoms)
    for h_group in hydrogen_groups:
        equivalent_groups.append(h_group['indices'])
    
    return equivalent_groups

#Write equivalent atoms to file in the required format
def write_equivalent_atoms(equivalent_groups, filename="similar.dat"):
    
    with open(filename, "w") as similar_file:
        for group in equivalent_groups:
            if len(group) > 1:
                # Use the first atom as reference
                reference = group[0]
                for atom in group[1:]:
                    similar_file.write(f"{atom+1:5d} {reference+1:5d}\n")  # Convert to 1-based indexing

def main(input_file, file_type):
    distance_file = "distance.dat"
    angle_file = "angle.dat"
    dihedral_file = "dihedral.dat"

    # Choose parsing method based on file_type
    is_fchk = file_type in [3, 4]
    if file_type == 2:
        coordinates, hessian, atomic_numbers, charges = parse_gau(input_file)
    else:
        coordinates, hessian, atomic_numbers, charges = parse_fchk(input_file)

    atom_pairs = read_distances(distance_file)
    angle_definitions = read_angles(angle_file)
    dihedral_definitions = read_dihedrals(dihedral_file)

    bonds = seminario_method(hessian, coordinates, atom_pairs, atomic_numbers)
    atom_names = [periodictable.elements[atomic_num].symbol for atomic_num in atomic_numbers]
    equivalent_groups = identify_equivalent_atoms_enhanced(coordinates, atomic_numbers, atom_names)
    write_equivalent_atoms(equivalent_groups)

    with open("bond_angle_dihedral_data.dat", "w") as file:
        # Write bond information
        file.write("\nBond Information:\n")
        for bond in bonds:
            atom1, atom2 = bond[0] - 1, bond[1] - 1
            calculated_distance = calculate_distance(coordinates[atom1], coordinates[atom2], is_fchk)
            force_constant = bond[2]
            file.write(f"{bond[0]:5d} {bond[1]:5d} {calculated_distance:10.3f} {force_constant:15.2f}\n")

        # Write angle information
        file.write("\nAngle Information:\n")
        for angle_def in angle_definitions:
            atom1, atom2, atom3, _ = angle_def
            idx1, idx2, idx3 = atom1 - 1, atom2 - 1, atom3 - 1
            calculated_angle = calculate_angle(coordinates[idx1], coordinates[idx2], coordinates[idx3])
            force_constant = calculate_angle_force_constant(hessian, coordinates, idx1, idx2, idx3, is_fchk)
            file.write(f"{atom1:5d} {atom2:5d} {atom3:5d} {calculated_angle:10.3f} {force_constant:15.2f}\n")

        # Write dihedral information
        file.write("\nDihedral Information:\n")
        for dihedral_def in dihedral_definitions:
            atom1, atom2, atom3, atom4, dihedral_value = dihedral_def
            file.write(f"{atom1:5d} {atom2:5d} {atom3:5d} {atom4:5d} {dihedral_value:10.3f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py input_file file_type")
        print("file_type: 2 for Gaussian output file, 3 or 4 for Gaussian fchk file")
        sys.exit(1)

    input_file = sys.argv[1]
    file_type = int(sys.argv[2])

    if file_type not in [2, 3, 4]:
        print("Error: file_type must be 2 (Gaussian output) or 3/4 (Gaussian fchk)")
        sys.exit(1)

    main(input_file, file_type)