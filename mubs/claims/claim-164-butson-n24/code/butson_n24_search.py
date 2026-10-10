import sys, os, time, json
import numpy as np

sys.path.insert(0, 'MUBs')
from butson_n24_vanishing_sum import cyclotomic_reductions
import mub_common as mc
import butson_exhaust as be

def ray_classes_numpy(n):
    reps_list = cyclotomic_reductions(n)
    reps = np.array(reps_list, dtype=np.int8)
    phi_n = reps.shape[1]
    slots_tail = 6 - 1
    
    Z = []
    M = []
    
    inner_grid = np.stack(np.meshgrid(*([np.arange(n)] * (slots_tail - 1)), indexing="ij"), axis=-1)
    inner_flat = inner_grid.reshape(-1, slots_tail - 1).astype(np.int16)
    inner_count = inner_flat.shape[0]

    for v1 in range(n):
        v1_col = np.full((inner_count, 1), v1, dtype=np.int16)
        z0_col = np.zeros((inner_count, 1), dtype=np.int16)
        vec = np.concatenate([z0_col, v1_col, inner_flat], axis=1)
        
        gathered = reps[vec]
        totals = gathered.sum(axis=1)
        z_mask = np.all(totals == 0, axis=1)
        
        for idx in np.nonzero(z_mask)[0]:
            Z.append(tuple(int(x) for x in vec[idx]))
            
        diffs = (vec[:, :, None] - vec[:, None, :]) % n
        diffs_reps = reps[diffs]
        norms = diffs_reps.sum(axis=(1, 2))
        m_hit = (norms[:, 0] == 6) & np.all(norms[:, 1:] == 0, axis=1)
        
        for idx in np.nonzero(m_hit)[0]:
            M.append(tuple(int(x) for x in vec[idx]))
            
    return Z, M

def exhaust_24():
    started = time.time()
    n = 24
    print("Computing ray classes...")
    Z, M = ray_classes_numpy(n)
    print(f"|Z| = {len(Z)}, |M| = {len(M)} (took {time.time()-started:.2f}s)")
    
    result = {
        "N": n,
        "zero_rays": len(Z),
        "mu_rays": len(M),
        "verdict": "no_quadruple",
    }
    
    zero_set = set(Z)
    mu_set = set(M)
    
    print("Building row graph and finding Hadamards...")
    row_graph = be.forward_graph(Z, zero_set, n)
    hadamard_row_sets = be.cliques(row_graph, 5)
    result["hadamards"] = len(hadamard_row_sets)
    print(f"Found {len(hadamard_row_sets)} Hadamards (took {time.time()-started:.2f}s)")
    
    census = {}
    total_completions = 0
    max_u = 0
    max_completions = 0
    zero_row = (0,) * 6
    
    for h_index, row_set in enumerate(hadamard_row_sets):
        h_rows = [zero_row] + [Z[i] for i in row_set]
        rays = be.completion_rays(h_rows, M, mu_set, n)
        max_u = max(max_u, len(rays))
        
        ray_graph = be.forward_graph(rays, zero_set, n)
        bases = be.cliques(ray_graph, 6)
        
        max_completions = max(max_completions, len(bases))
        total_completions += len(bases)
        
        key = "%d/%d" % (len(rays), len(bases))
        census[key] = census.get(key, 0) + 1
        
        for a in range(len(bases)):
            for b in range(a + 1, len(bases)):
                if be.compatible_completions(bases[a], bases[b], rays, mu_set, n):
                    result["verdict"] = "quadruple"
                    result["witness"] = be.exact_hit_check(n, h_rows, rays, bases[a], bases[b])
                    result["seconds"] = round(time.time() - started, 3)
                    print(f"Found quadruple!")
                    return result
        
        if (h_index + 1) % 100 == 0:
            print(f"Checked {h_index + 1}/{len(hadamard_row_sets)} Hadamards...")
                    
    result.update({
        "total_completions": total_completions,
        "max_u": max_u,
        "max_completions": max_completions,
        "types": dict(sorted(census.items())),
        "reason": "all_completion_pairs_excluded",
        "seconds": round(time.time() - started, 3)
    })
    
    print(f"Finished in {time.time()-started:.2f}s")
    return result

def verify_n24_log():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_path = os.path.join(script_dir, "../logs/butson_n24_numpy_search.json")
    assert os.path.exists(log_path), f"Log file not found at {log_path}"
    with open(log_path) as f:
        committed = json.load(f)
    print("Re-running exhaust_24() for independent verification (expect ~4.5 min)...")
    fresh = exhaust_24()
    # The 'seconds' timing field is intrinsically non-deterministic; every
    # other field must match the committed log byte for byte under the
    # canonical sort_keys / separators serialization.
    for d in (committed, fresh):
        d.pop("seconds", None)
    def canon(d):
        return json.dumps(d, sort_keys=True, separators=(",", ":"))
    if canon(fresh) != canon(committed):
        raise AssertionError(
            "committed search log differs from fresh exhaust_24() output:\n"
            f"  committed: {canon(committed)}\n"
            f"  fresh:     {canon(fresh)}"
        )
    print("VERIFIED: committed Butson N=24 search log reproduces exhaust_24() exactly.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Butson N=24 MUB Quadruple Exhaustion Search")
    parser.add_argument("--verify", action="store_true", help="Verify the committed N=24 search log")
    args = parser.parse_args()
    
    if args.verify:
        verify_n24_log()
    else:
        res = exhaust_24()
        log_path = "problems/mub-quad/logs/butson_n24_numpy_search.json"
        with open(log_path, "w") as f:
            json.dump(res, f, indent=2, sort_keys=True)

