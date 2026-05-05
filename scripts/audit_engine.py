import struct
import math
import subprocess
import re
from pathlib import Path

def read_csr(filename):
    if not Path(filename).is_file():
        raise FileNotFoundError(f"Missing CSR export: {filename}")
    with open(filename, "rb") as f:
        fmt_header = "QQ"
        header_size = struct.calcsize(fmt_header)
        header = f.read(header_size)
        if len(header) != header_size:
            raise ValueError(f"CSR file {filename} is too short to contain a header")
        n_states, n_transitions = struct.unpack(fmt_header, header)
        row_bytes = f.read(8 * (n_states + 2))
        out_bytes = f.read(8 * (n_transitions + 1))
        edge_bytes = f.read(8 * (n_transitions + 1))
        if len(row_bytes) != 8 * (n_states + 2):
            raise ValueError(f"CSR file {filename} has a truncated row pointer block")
        if len(out_bytes) != 8 * (n_transitions + 1):
            raise ValueError(f"CSR file {filename} has a truncated out-state block")
        if len(edge_bytes) != 8 * (n_transitions + 1):
            raise ValueError(f"CSR file {filename} has a truncated edge block")
        row_ptr = list(struct.unpack("Q" * (n_states + 2), row_bytes))
        out_states = list(struct.unpack("Q" * (n_transitions + 1), out_bytes))
        edges = list(struct.unpack("Q" * (n_transitions + 1), edge_bytes))
    return n_states, row_ptr, out_states, edges

def read_vector(filename):
    if not Path(filename).is_file():
        raise FileNotFoundError(f"Missing eigenvector export: {filename}")
    with open(filename, "r") as f:
        return [float(line.strip()) for line in f]

def parse_metric(output, pattern, label):
    match = re.search(pattern, output)
    if not match:
        raise RuntimeError(f"Could not parse {label} from tm_master output")
    return float(match.group(1))

def matvec(x, v, n_states, row_ptr, out_states, edges):
    res = [0.0] * n_states
    for i in range(1, n_states + 1):
        start = row_ptr[i]
        end = row_ptr[i+1]
        for idx in range(start + 1, end + 1):
            out = out_states[idx]
            w = edges[idx]
            res[i-1] += (x**w) * v[out-1]
    return res

def get_rho(x, n_states, row_ptr, out_states, edges):
    v = [1.0] * n_states
    for _ in range(200):
        v_new = matvec(x, v, n_states, row_ptr, out_states, edges)
        rho = max(v_new)
        if rho == 0: return 0
        v = [vi/rho for vi in v_new]
    return rho

def audit_full(L, M, mode, x0_code, alpha_code, beta_code):
    mode_suffixes = ["std", "ham", "2sap", "2sap_ham"]
    suffix = mode_suffixes[mode] if 0 <= mode < len(mode_suffixes) else "unk"
    
    csr_file = f"data/TransferMatrix/TMresults/CSR_L{L}M{M}_{suffix}.bin"
    lv_file = f"data/TransferMatrix/TMresults/L_Evector_L{L}M{M}_{suffix}.txt"
    rv_file = f"data/TransferMatrix/TMresults/R_Evector_L{L}M{M}_{suffix}.txt"
    
    n_states, row_ptr, out_states, edges = read_csr(csr_file)
    eta = read_vector(lv_file) 
    xi = read_vector(rv_file)  
    
    LmultR = sum(e * x for e, x in zip(eta, xi))
    g_xi = matvec(x0_code, xi, n_states, row_ptr, out_states, edges)
    res_inf = max(abs(a - b) for a, b in zip(g_xi, xi))
    
    B = 0
    for i in range(1, n_states + 1):
        start = row_ptr[i]
        end = row_ptr[i+1]
        for idx in range(start + 1, end + 1):
            out = out_states[idx]
            w = edges[idx]
            B += eta[i-1] * w * (x0_code**w) * xi[out-1]
    
    s_analytical = B / (LmultR * x0_code)
    beta_spec = x0_code * s_analytical
    alpha_spec = 1.0 / s_analytical
    
    print(f"\n--- AUTHORITATIVE AUDIT REPORT (LxM={L}x{M}) ---")
    print(f"x0: {x0_code:.12f}")
    print(f"Spectral Residual: {res_inf:.2e}")
    print(f"s (Analytical):    {s_analytical:.12f}")
    
    print(f"\n--- Mismatch Analysis ---")
    print(f"{'Value':<10} | {'Code Output':<15} | {'Spec Target':<15} | {'Rel Error':<10}")
    print("-" * 60)
    print(f"{'Beta':<10} | {beta_code:15.6e} | {beta_spec:15.6e} | {abs(beta_code - beta_spec)/beta_spec:10.2e}")
    print(f"{'Alpha':<10} | {alpha_code:15.6e} | {alpha_spec:15.6e} | {abs(alpha_code - alpha_spec)/alpha_spec:10.2e}")
    
    print(f"\n--- Pole Residue Convergence (phi-phi) ---")
    target_res = alpha_spec * (xi[0] * eta[0] / LmultR)
    print(f"Theoretical Residue: {target_res:.8e}")
    
    for eps in [1e-3, 1e-4, 1e-5]:
        x = x0_code - eps
        inv_e1_1 = 0
        v = [0.0] * n_states
        v[0] = 1.0
        for k in range(5000): # More iterations for eps=1e-5
            inv_e1_1 += v[0]
            v = matvec(x, v, n_states, row_ptr, out_states, edges)
            if max(v) < 1e-15: break
        actual_res = eps * inv_e1_1
        print(f"eps={eps:.0e}: Res={actual_res:.8e} (Err={abs(actual_res - target_res)/target_res:.2e})")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Audit TMcalc Research Results")
    parser.add_argument("-L", type=int, default=2, help="Lattice Width")
    parser.add_argument("-M", type=int, default=1, help="Lattice Height")
    parser.add_argument("-m", "--mode", type=int, default=0, help="Simulation Mode (0:Std, 1:Ham, 2:2SAP, 3:2SAP_Ham)")
    
    args = parser.parse_args()
    L, M, mode = args.L, args.M, args.mode
    
    cmd = ["bin/tm_master", "-L", str(L), "-M", str(M), "-m", str(mode)]
    print(f"Running Engine: {' '.join(cmd)}")
    output = subprocess.check_output(cmd, text=True)
    
    x0 = parse_metric(output, r"x_0\):\s+([\d\.]+)", "critical fugacity")
    alpha = parse_metric(output, r"Alpha \(Amplitude\):\s+([\d\.e\-\+]+)", "alpha")
    beta = parse_metric(output, r"Beta \(Growth Parameter\):\s+([\d\.e\-\+]+)", "beta")
    
    audit_full(L, M, mode, x0, alpha, beta)
